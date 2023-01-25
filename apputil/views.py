from django.contrib.admin.views.decorators import staff_member_required
from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView, CreateView, DeleteView

from .forms import ApplicationUser_form, Dictionary_form, Login_form
from .models import ApplicationUser, Dictionary
from .utils_dataimport import import_excel
from dorganism.models import Organism, Taxonomy
from apputil.utils import FilteredListView
from dorganism.utils import Dictionaryfilter
from adjcoadd.constants import *
from ddrug.models import VITEK_Card, VITEK_ID, VITEK_AST

# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")

## =================================APP Home========================================

# import setup
@login_required(login_url='/')
def index(req):
    # print(setup.version)
    object_1=Organism.objects.count()
    object_2=Taxonomy.objects.count()
    return render(req, 'dorganism/home.html', {'objects_org': object_1, 'objects_taxo':object_2,})
## =================================APP Home======================================##

## =================================APP Log in/out =================================
def login_user(req):
    if req.user.is_authenticated:
        return redirect("index")
    else:

        if req.method=='POST':
            form=Login_form(data=req.POST)
            username_ldap=req.POST.get('username')
        # print(user)
            if form.is_valid():
                print("form is valid")
                user=form.get_user()
                login(req, user, backend="django_auth_ldap.backend.LDAPBackend",)
                return redirect("index")
        
            else:
                messages.warning(req, ' no permission for this application, please contact Admin!')
                return redirect("/")
        else:
            form = Login_form()
        return render(req, 'registration/login.html', {'form': form})    

def logout_user(req):
    logout(req)    
    return redirect("/")

# =================================APP Log in/out ==================================##

## =========================Application Users View====================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = '/'

    def test_func(self):
        return self.request.user.has_permission('Admin')

@login_required(login_url='/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    return render(req, 'apputil/userprofile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, ListView):
    login_url = '/'
    model=ApplicationUser
    fields='__all__'
    template_name = 'apputil/appUsers.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context
        
from django.http import QueryDict
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def updateApplicationuser(req, pk):
    object_=get_object_or_404(ApplicationUser, pk=pk)
    form=ApplicationUser_form(instance=object_)
    context={
        "form":form,
        "object":object_,
        }
    if req.method=='PUT':
        qd=QueryDict(req.body).dict()       
        form=ApplicationUser_form(data=qd, instance=object_)
        if form.is_valid():
            instance=form.save()
            context={"object":object_}
            return render(req, "apputil/appuser_tr.html", context)
    return render(req, "apputil/appUsersUpdate.html", context)

class AppUserCreateView(SuperUserRequiredMixin, CreateView):
    model=ApplicationUser
    fields=['name', 'username', 'permission',]
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')


class AppUserDeleteView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    template_name='apputil/appUsersDel.html'
    success_url = reverse_lazy('userslist')
    fields=['is_appuser']

    def form_valid(self, form):

        form.instance.is_appuser==False
        print(form.instance.is_appuser)
        return super().form_valid(form)
# =========================Application Users View==================================##

## ========================Dictionary View===========================================

class DictionaryView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Dictionary
    template_name='apputil/dictList.html'
    filterset_class = Dictionaryfilter
    model_fields=DICTIONARY_FIELDs


# 
@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def createDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    form=Dictionary_form()
    form_error=False
    if req.method=='POST':
        form=Dictionary_form(req.POST)
        try:
            if form.is_valid:
                print("form is valid")
                instance=form.save(commit=False)
                instance.save(**kwargs)
                return redirect("dict_view")
        except Exception as err:
       
            form_error=True    
            messages.error(req, form.errors)
            return redirect("dict_view")
    
    return render(req, 'apputil/dictCreate.html', {'form': form, 'form_error':form_error})
## ============================Dictionary View======================================##
@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def updateDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
        dict_value=req.POST.get("dict_value")
        dict_class=req.POST.get("dict_class")
        dict_desc=req.POST.get("dict_desc")
        object_=get_object_or_404(Dictionary, dict_value=dict_value)
        try:
            if object_:
                object_.dict_class=dict_class
                object_.dict_desc=dict_desc
                print(object_.dict_desc)
                object_.save(**kwargs)
                return JsonResponse({"result": "Saved"})
        except Exception as err:
            return JsonResponse({"result": err})
    
    return JsonResponse({})


# ==========================File Process=====================================

from django.core.files.storage import FileSystemStorage
from django.views import View
from django import forms
import json
from django.core import serializers
import os
from .utils_dataimport import FileValidator, uploadedfile_process
from django.core.exceptions import ValidationError
from django.db import transaction, IntegrityError
from pathlib import Path
from django.conf import settings
from .utils import instance_dict, Validation_Log
from asgiref.sync import sync_to_async


# set filefield Validator
# validate_file = FileValidator(#max_size=1024 * 100, 
#                              content_types=('text/csv', 'application/pdf','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
# create array for files if infected
# infected_files = []
# setup unix socket to scan stream
# cd = clamd.ClamdUnixSocket()

if settings.DEVELOPMENT:
    path='uploads'
else:
    Base_dir = Path(__file__).resolve().parent.parent.parent
    path=os.path.abspath(os.path.join(Base_dir, 'uploads'))

    # #delete task

def delete_file(file_path):
    file_name=file_path.split("/")[2]
    print(file_name)
    file_full_path=os.path.join(settings.MEDIA_ROOT, file_name)
    print(file_full_path)
    try:
        os.unlink(file_full_path)
        print("removed!")
    except Exception as err:
        print(err)


class FileUploadForm(SuperUserRequiredMixin, forms.Form):
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True,}), )#validators=[validate_file])
    

class Importhandler_VITEK(SuperUserRequiredMixin, View):
    
    form_class=FileUploadForm
    file_url=[]
    data_list=[]
    data_model='default'
 
    def get(self, request):
        form = self.form_class
        for f in os.listdir(path):
            print(f)
        return render(request, 'ddrug/importdata_vitek.html', { 'form': form, })
    
   

    
    def post(self, request):
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form
        vLog = Validation_Log('VITEK PDF')
        kwargs={}
        kwargs['user']=request.user
        vCards=[]
        vID=[]
        vAST=[]
        self.data_model=request.POST.get('file_data')
        myfiles=request.FILES.getlist('file_field')
        self.file_url=[]
        try:
            if form.is_valid():
                print(myfiles)
                # scan_results = cd.instream(myfile) # scan_results['stream'][0] == 'OK' or 'FOUND'
                for f in myfiles:
                    fs=FileSystemStorage()
                    filename=fs.save(f.name, f)
                    self.file_url.append(fs.url(filename))
                    print(self.file_url)
                    try:
                        vCards,vID,vAst=import_excel(fs.url(filename), self.data_model)
                        print("file checked")
           
                    except Exception as err:
                        print(f"uploaderror is {err}")

                    # return messages.warning(request, f'There is {err} error, upload again')
                context['file_pathlist']=self.file_url
                context['data_model']=self.data_model
                return render(request,'ddrug/importdata_vitek.html', context)
            else:
                messages.warning(request, f'There is {form.errors} error, upload again')          

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')
        
        if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
            print(f'self.file_url is : {self.file_url}')
            table_name=[","]
            validate_result=[","]
            file_report=[","]
            process_name=request.POST.get('type')
            file_pathlist=request.POST.getlist("filepathlist[]")
            print(f'selected : {file_pathlist}')
            data_model=request.POST.get("datamodel")
                # uploadedfile_process(request, table_name,validate_result, file_report, process_name, f, data_model, vLog)
            if process_name=='Validation':
                for f in file_pathlist:
                # models objects list coming from parsed file
                    vCards,vID,vAst=import_excel(f, data_model)
                # validating each objectslist 
                    if vCards:
                        table_name.append("Vitek_card")
                        for e in vCards:
                            djCard=VITEK_Card.check_from_dict(e, vLog)
                            validate_result.append(f"CARD status: {djCard.validStatus}")
                            file_report.append(str(vLog.show()))
        
                    if vID:
                        table_name.append("Vitek_id")
                        for e in vID:
                            djID=VITEK_ID.check_from_dict(e, vLog)
                            validate_result.append(f"ID status: {djID.validStatus}")
                            file_report.append(str(vLog.show()))
        
                    if vAst:
                        table_name.append("Vitek_ast")
                        for e in vAst:
                            djAst=VITEK_AST.check_from_dict(e, vLog)
                            validate_result.append(f"AST status: {djAst.validStatus}")
                            file_report.append(str(vLog.show()))               
        
                return JsonResponse({"table name":"VITEK".join(table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report)})                                   
       
            elif process_name=='Cancel':
                # Cancel Task
                for f in file_pathlist:
                    delete_file(file_path=f)
                return JsonResponse({"table name":(",").join( table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report)})                                   
    
            elif process_name=='DB_Validation':
                # import data to DB
                for f in file_pathlist:             
                    vCards, vID, vAst=import_excel(f, data_model)
        # validating each objectslist 
                    if vCards:
                        table_name.append("Vitek_card")
                        for e in vCards:
                            djCard=VITEK_Card.check_from_dict(e, vLog)
                            if djCard.validStatus:
                                try:
                                    djCard.save(**kwargs)
                                except Exception as err:
                                    validate_result.append(f"catch Exception CARD {err}")
                            validate_result.append(f"CARD status: {djCard.validStatus}")
                            file_report.append(str(vLog.show()))
                
                    if vID:
                        table_name.append("Vitek_id")
                        for e in vID:
                            djID=VITEK_ID.check_from_dict(e, vLog)
                            if djID.validStatus:
                                try:
                                    djID.save(**kwargs)
                                except Exception as err:
                                    validate_result.append(f"catch Exception ID {err}") 
                            validate_result.append(f"ID status: {djID.validStatus}")
                            file_report.append(str(vLog.show()))
                
                    if vAst:
                        table_name.append("Vitek_ast")
                        for e in vAst:
                            djAst=VITEK_AST.check_from_dict(e, vLog)
                            if djAst.validStatus:
                                try:
                                    djAst.save(**kwargs)
                                except Exception as err:
                                    validate_result.append(f"catch Exception Ast {err}")    
                            validate_result.append(f"AST status: {djAst.validStatus}")
                            file_report.append(str(vLog.show()))               
                
                return JsonResponse({"table name":"VITEK".join(table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report), 'status':"Data Saved!"})

           
        return render(request, 'ddrug/importdata_vitek.html', context)

async_function = sync_to_async(Importhandler_VITEK.get, thread_sensitive=False)
async_function = sync_to_async(Importhandler_VITEK.post, thread_sensitive=False)

      


  