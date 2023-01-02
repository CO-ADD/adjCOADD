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
from adjcoadd.utils import import_excel
from dorganism.models import Organism, Taxonomy

# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")



## =================================APP Home========================================
@login_required(login_url='/')
def index(req):
    object_1=Organism.objects.count()
    object_2=Taxonomy.objects.count()
    return render(req, 'dorganism/home.html', {'objects_org': object_1, 'objects_taxo':object_2})
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
        return self.request.user.is_superuser

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

class AppUserCreateView(SuperUserRequiredMixin, CreateView):
    model=ApplicationUser
    fields=['name', 'username', 'permission']
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context

class AppUserUpdateView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    fields=['name', 'permission', ]
    template_name = 'apputil/appUsersUpdate.html'
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

class DictionaryView(LoginRequiredMixin, ListView):
    login_url = '/'
    model=Dictionary
    fields="__all__"
    template_name='apputil/dictList.html'
    
    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.filter(astatus__gte=0)
        return context
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



# ==========================File Process=====================================

from django.core.files.storage import FileSystemStorage
from django.views import View
from django import forms
import json
from django.core import serializers
import os
from .forms import FileValidator
from django.core.exceptions import ValidationError
from django.db import transaction, IntegrityError

# set filefield Validator
# validate_file = FileValidator(max_size=1024 * 100, 
#                              content_types=('text/csv', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
# create array for files if infected
# infected_files = []
# setup unix socket to scan stream
# cd = clamd.ClamdUnixSocket()

class FileUploadForm(forms.Form):
    file_data= forms.ChoiceField(choices=(('Taxonomy', 'Taxonomy'),('Organism', 'Organism'),))
    file_field = forms.FileField()#(validators=[validate_file])

class Importhandler(View):
    template_name='apputil/importdata.html'
    form_class=FileUploadForm
    file_url=''
    data_list=[]
    data_model='default'

    def get(self, request):
        form = self.form_class
        return render(request, 'apputil/importdata.html', { 'form': form, })

    def post(self, request):
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form

        try:
            if form.is_valid():
                myfile=request.FILES['file_field']
                self.data_model=request.POST.get('file_data')
                # scan_results = cd.instream(myfile) # scan_results['stream'][0] == 'OK' or 'FOUND'
                fs=FileSystemStorage()
                filename=fs.save(myfile.name, myfile)
                self.file_url=fs.url(filename)
                context['file_path']=self.file_url
                context['data_model']=self.data_model
                return render(request,'apputil/importdata.html', context)
            else:
                messages.warning(request, f'There is {form.errors} error, upload again')          

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')

        
        if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
            process_name=request.POST.get('type')
            print(f'step is {process_name}')
            #type: Validation, Cancel, DB_Validation, Confirm, RollBack, Save-Data
            if process_name=='Validation':
                file_path=request.POST.get("filepath")
                data_model=request.POST.get("datamodel")
                Importhandler.data_list=import_excel(file_path, data_model)           
                if isinstance(Importhandler.data_list, Exception):
                    result=str(Importhandler.data_list)
                    status='Form Errors'
                else:
                    result='pass'
                    status='Form is Valid'
                return JsonResponse({"task_user": request.user.pk, 'task_status':status, 'task_result': result})
            elif process_name=='Cancel':
                # Cancel Task
                Importhandler.delete_task(request)
                return JsonResponse({})
            elif process_name=='DB_Validation':
                # import data to DB
                print(f'save data {Importhandler.data_list} to db')
                errList=[]
                #  save to db function
                
                for obj in Importhandler.data_list:
                    print(obj.pk)
                    if Taxonomy.objects.filter(pk=obj.pk):
                        return JsonResponse({'status':'DATA exists'}) 
                    try:
                        obj.save(commit=False)
                        print('save')
                    except Exception as err:
                        print(err)
                        errList.append[err]
                        return JsonResponse({"status":errList})
            
                return JsonResponse({"status":"SUCCESS"}, status=200)
                  
            elif process_name=='Save-Data':
                print(Importhandler.data_list)
                with transaction.atomic(using='dorganism'):
                    for obj in Importhandler.data_list:
                        try:
                            obj.save()
                            print('save')
                        except Exception as err:
                            print(err)
                            errList.append[err]
                            return JsonResponse({"status":errList})
                return JsonResponse({"status":"Data Saved!"}, status=200)
        return render(request, 'apputil/importdata.html', context)
    
    #delete task
    @csrf_exempt
    def delete_task(request):
        print(Importhandler.file_url)
        for filename in os.listdir("uploads"):
                file_path=os.path.join("uploads", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)
        return JsonResponse({})


  