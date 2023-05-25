import os
from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse, QueryDict
from django.http import JsonResponse, QueryDict
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView, CreateView, DeleteView

from adjcoadd.constants import *
from dorganism.models import Organism, Taxonomy
from ddrug.models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub
from dgene.models import Gene, WGS_CheckM, WGS_FastQC, ID_Pub, ID_Sequence

from apputil.forms import AppUserfilter, Dictionaryfilter, ApplicationUser_form, Dictionary_form, Login_form
from apputil.models import ApplicationUser, Dictionary
from apputil.utils.views_base import SuperUserRequiredMixin, permission_not_granted, SimplecreateView, HtmxupdateView
from apputil.utils.filters_base import FilteredListView
from apputil.utils.files_upload import Importhandler

## =================================APP Home========================================

# import setup
@login_required(login_url='/')
def index(req):

    nDict = {    
        'nOrg':    Organism.objects.count(),
        'nTax':    Taxonomy.objects.count(),
        'nDrug':   Drug.objects.count(),
        'nVCard':  VITEK_Card.objects.count(),
        'nVID':    VITEK_ID.objects.count(),
        'nVAST':   VITEK_AST.objects.count(),
        'nMICC':   MIC_COADD.objects.count(),
        'nMICP':   MIC_Pub.objects.count(),
#        'nBP':     Breakpoints.objects.count(),
        'nGene':   Gene.objects.count(),
        'nCheckM': WGS_CheckM.objects.count(),
        'nFastQC': WGS_FastQC.objects.count(),
        'nIDP':    ID_Pub.objects.count(),
        'nIDS':    ID_Sequence.objects.count(),
    }
    return render(req, 'home.html', nDict)

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

@login_required(login_url='/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    return render(req, 'apputil/appUserProfile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=ApplicationUser
    template_name = 'apputil/appUsers.html'  
    filterset_class = AppUserfilter
    model_fields=model.HEADER_FIELDS

class AppUserCreateView(SuperUserRequiredMixin, SimplecreateView):
    form_class = ApplicationUser_form
    template_name = 'apputil/appUsersCreate.html'

    def post(self, request, *args, **kwargs):
        form =self.form_class(request.POST)
        if form.is_valid():
            instance=form.save()
            return redirect(request.META['HTTP_REFERER'])
        else:
            messages.error(request, form.errors)
            return redirect(request.META['HTTP_REFERER'])

##
## here used HTMX
from apputil.utils.views_base import HtmxupdateView
class ApplicationUserUpdateView(HtmxupdateView):
    form_class=ApplicationUser_form
    template_name="apputil/appUsersUpdate.html"
    template_partial="apputil/appuser_tr.html"
    model=ApplicationUser

    def put(self, request, *args, **kwargs):
        pk=kwargs.get("pk")
        object_=self.get_object(pk)
        qd=QueryDict(request.body).dict()
        form =self.form_class(data=qd, instance=object_)
        context={
        "form":form,
        "object":object_,
    }
        if form.is_valid():           
            object_new=form.save()                
            return render(request, self.template_partial, context)
        else:
            messages.error(request, form.errors)
            return render(request, self.template_partial, context)



class AppUserDeleteView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    template_name='apputil/appUsersDel.html'
    success_url = reverse_lazy('userslist')
    fields=['is_appuser']

    def form_valid(self, form):
        form.instance.is_appuser==False
        return super().form_valid(form)
# =========================Application Users View==================================##

## ========================Dictionary View===========================================

class DictionaryView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Dictionary
    template_name='apputil/dictList.html'
    filterset_class = Dictionaryfilter
    model_fields=model.HEADER_FIELDS

    
# 
class DictionaryCreateView(SuperUserRequiredMixin, SimplecreateView):
    form_class=Dictionary_form
    template_name='apputil/dictCreate.html'

## ============================Dictionary View======================================##
@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def updateDictionary(req):
    kwargs={'user': req.user}
    if req.user.has_permission('Admin'):
        if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
            dict_value=req.POST.get("dict_value") 
            type=req.POST.get("type") or None
            value=req.POST.get("value").strip() or None
            object_=get_object_or_404(Dictionary, dict_value=dict_value)
            try:
                if object_:
                    if type=='dict_class':
                        object_.dict_class=value
                    if type=='dict_desc':
                        object_.dict_desc=value
                    object_.save(**kwargs)
                    return JsonResponse({"result": "Saved"})
            except Exception as err:
                return JsonResponse({"result": err})
    else:
        return JsonResponse({"result": "Permission_denied"})
    return JsonResponse({})


@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def deleteDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
        dict_value=req.POST.get("dict_value") 
        object_=get_object_or_404(Dictionary, dict_value=dict_value)
        print(object_)
        try:
            if object_:
                object_.delete(**kwargs)
                return JsonResponse({"success": 'data deleted'})
        except Exception as err:
            return JsonResponse({"error": err})
    
    return JsonResponse({})


############################################### Export CSV View ###########################################
import csv
import datetime 
from django.apps import apps

@login_required
@user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
def exportCSV(request):
    if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
        data_arr = request.POST.getlist('data_arr[]')
        data_fields = request.POST.getlist('fields[]')
        model_name=request.POST.get('model_name')

        try:
            model=apps.get_model('dorganism', model_name)
        except:
            model=apps.get_model('ddrug', model_name)
        query=model.objects.filter(pk__in=data_arr)
        response = HttpResponse(content_type='text/csv')
        file_name = "fltred_loaction_data" + str(datetime.date.today()) + ".csv"
        writer = csv.writer(response)
        writer.writerow(data_fields)
        for i in query.values_list(*data_fields):
            writer.writerow(i)
        response['Content-Disposition'] = 'attachment; filename = "' + file_name + '"'
        return response


# =============Import Dictionary and appUsers via Excel==================
from .utils.files_upload import FileUploadForm, OverwriteStorage, file_location
from .utils.validation_log import Validation_Log
# from impdata.a_upload_AppUtil import update_AppUser_xls, update_Dictionary_xls
# from impdata.c_upload_dDrug import update_Drug_xls

class Importhandler_apputils(Importhandler):
    form_class=FileUploadForm
    template_name='apputil/importhandler_excel.html'
    upload_model_type=None
    lObject={}   # store results parsed by all uploaded excel files with key-filename, value-parsed result array
    # vLog = Validation_Log("Vitek-pdf")
    
    def post(self, request, process_name):
        process_name=process_name
        uploadDir=file_location(request) # define file store path during file process
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form
        context["process_name"]=process_name
        
       
        kwargs={}
        kwargs['user']=request.user
        
        self.file_url=[]
        self.data_model=request.POST.get('file_data') or None
        myfiles=[request.FILES.get('file_field'),]
        print(myfiles)
      
        try:
        # Uploading Verifying     
            if form.is_valid():
                self.lObject.clear()

        # Uploading Parsing 
                for f in myfiles:
                    fs=OverwriteStorage(location=uploadDir)
                    xlFile=fs.save(f.name, f)
                    uploadFile = os.path.join(uploadDir, xlFile)
                   
                    try:
                        if process_name=="Dictioanry":
                            # update_AppUser_xls(uploadFile, XlsSheet="Dictionary", upload=True, uploaduser=request.user, lower=False)
                            context["excel_upload_info"]="Saved in Dictionary Datatable"
                        elif process_name=="ApplicationUser":
                            # update_Dictionary_xls(uploadFile,XlsSheet="User",upload=True)
                            context["excel_upload_info"]="Saved in AppUser Datatable"
                        elif process_name=="Drug":
                            # update_Drug_xls(uploadFile, XlsSheet="Drug", upload=False, uploaduser=request.user, lower=True)
                            context["excel_upload_info"]="Saved in Drug Datatable"
                            
                    except Exception as err:
                        context["excel_upload_info"]=f"{err} cause failed import"
                        self.delete_file(file_name=filename)
                        return render(request, template_name, context)
                  

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')
           
        return render(request, self.template_name, context)
