import os
from formtools.wizard.views import SessionWizardView
from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse, QueryDict
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView
from django.views.generic.detail import DetailView
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from dorganism.models import Organism, Taxonomy
from ddrug.models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
from dgene.models import Gene, WGS_CheckM, WGS_FastQC, ID_Pub, ID_Sequence

from apputil.forms import AppUserfilter, Dictionaryfilter, ApplicationUser_form, Dictionary_form, Login_form, Image_form, Document_form 
from apputil.models import ApplicationUser, Dictionary, Image, Document
from apputil.utils.views_base import SuperUserRequiredMixin, permission_not_granted, SimplecreateView, SimpleupdateView,SimpledeleteView, HtmxupdateView, CreateFileView
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
        'nBP':     Breakpoint.objects.count(),
        'nGene':   Gene.objects.count(),
        'nCheckM': WGS_CheckM.objects.count(),
        'nFastQC': WGS_FastQC.objects.count(),
        'nIDP':    ID_Pub.objects.count(),
        'nIDS':    ID_Sequence.objects.count(),
    }
    return render(req, 'home.html', nDict)

# 
# Handler 404
def custom_page_not_found_view(request, exception):
   
    return render(request, 'registration/error404.html')
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

## =========================Application Users View====================================

class AppUserDetailView(DetailView):
    model = ApplicationUser
    template_name = 'apputil/appUserProfile.html' 

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # add extra context here...
        return context

@login_required(login_url='/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    return render(req, 'apputil/appUserProfile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = ApplicationUser
    template_name = 'apputil/appUsers.html'  
    filterset_class = AppUserfilter
    model_fields=model.HEADER_FIELDS

class AppUserCreateView(SuperUserRequiredMixin, SimplecreateView):
    form_class = ApplicationUser_form
    template_name = 'apputil/appUsersCreate.html'

    def post(self, request, *args, **kwargs):
        form = self.form_class(request.POST)
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
    form_class = ApplicationUser_form
    template_name = "apputil/appUsersUpdate.html"
    template_partial = "apputil/appuser_tr.html"
    model = ApplicationUser

    def put(self, request, *args, **kwargs):
        pk = kwargs.get("pk")
        object_= self.get_object(pk)
        qd = QueryDict(request.body).dict()
        form = self.form_class(data=qd, instance=object_)
        context = {
        "form":form,
        "object":object_,
    }
        if form.is_valid():           
            form.save()                
            return render(request, self.template_partial, context)
        else:
            messages.error(request, form.errors)
            return render(request, self.template_partial, context)

class AppUserDeleteView(SuperUserRequiredMixin, UpdateView):
    model = ApplicationUser
    template_name = 'apputil/appUsersDel.html'
    success_url = reverse_lazy('userslist')
    fields=['is_appuser']

    def form_valid(self, form):
        form.instance.is_appuser==False
        return super().form_valid(form)




## ========================Dictionary View===========================================

class DictionaryView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = Dictionary
    template_name = 'apputil/dictList.html'
    filterset_class = Dictionaryfilter
    model_fields = model.HEADER_FIELDS

    
# 
class DictionaryCreateView(SuperUserRequiredMixin, SimplecreateView):
    form_class = Dictionary_form
    template_name = 'apputil/dictCreate.html'


@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def updateDictionary(req):
    kwargs = {'user': req.user}
    if req.user.has_permission('Admin'):
        if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
            dict_value = req.POST.get("dict_value") 
            type = req.POST.get("type") or None
            value = req.POST.get("value").strip() or None
            object_= get_object_or_404(Dictionary, dict_value=dict_value)
            try:
                if object_:
                    if type == 'dict_class':
                        object_.dict_class=value
                    if type == 'dict_desc':
                        object_.dict_desc=value
                    if type == 'dict_sort':
                        object_.dict_sort=int(value)
                    object_.save(**kwargs)
                    return JsonResponse({"result": "Saved"})
            except Exception as err:
                return JsonResponse({"result": err})
    else:
        return JsonResponse({"result": "Permission_denied"})
    return JsonResponse({})


@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def deleteDictionary(req):
    kwargs = {}
    kwargs['user'] = req.user
    if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
        dict_value = req.POST.get("dict_value") 
        object_ = get_object_or_404(Dictionary, dict_value=dict_value)
        try:
            if object_:
                object_.delete(**kwargs)
                return JsonResponse({"success": 'data deleted'})
        except Exception as err:
            return JsonResponse({"error": err})
    
    return JsonResponse({})

## ============================File and Image View======================================## 

class CreatedocumentView(CreateFileView):
    form_class = Document_form
    model = Organism
    file_field = 'doc_file'
    related_name = 'assoc_documents'
    transaction_use_manytomany = 'dorganism'

class CreateimageView(CreateFileView):
    form_class = Image_form
    model = Organism
    file_field = 'image_file'
    related_name = 'assoc_images'
    transaction_use_manytomany = 'dorganism'
    
class DocDeleteView(SimpledeleteView):
    model = Document

class ImageDeleteView(SimpledeleteView):
    model = Image
# =========================== Export CSV View =============================
from .utils.views_base import DataExportBaseView
class DataExportView(DataExportBaseView):
    pass

# =============Import Dictionary and appUsers via Excel==================
from .utils.files_upload import FileUploadForm, OverwriteStorage, file_location
from .utils.validation_log import Validation_Log
from .utils.data_style import convert_heatmap
from .utils.form_wizard_tools import SelectFile_StepForm

class Importhandler_apputils(Importhandler):
    form_class=SelectFile_StepForm
    template_name='apputil/importhandler_excel.html'
    upload_model_type=None
    lObject={}   
    
    def post(self, request, process_name):
        process_name=process_name
        uploadDir=file_location(instance=request.user) # define file store path during file process
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form
        context["process_name"]=process_name  
       
        kwargs={}
        kwargs['user']=request.user
        myfiles=[]
        self.file_url=[]
        self.data_model=request.POST.get('file_data') or None
        myfiles.extend(request.FILES.getlist('multi_files'),)
        print(myfiles)
      
        try:
            for f in myfiles:
                fs=OverwriteStorage(location=uploadDir)
                xlFile=fs.save(f.name, f)
                uploadFile = os.path.join(uploadDir, xlFile)
                print(f'uploadfile: {uploadFile}')
               
                # try:
                if process_name=="Dictioanry":
                    # update_AppUser_xls(uploadFile, XlsSheet="Dictionary", upload=True, uploaduser=request.user, lower=False)
                    context["excel_upload_info"]="Saved in Dictionary Datatable"
                elif process_name=="ApplicationUser":
                    # update_Dictionary_xls(uploadFile,XlsSheet="User",upload=True)
                    context["excel_upload_info"]="Saved in AppUser Datatable"
                elif process_name=="Drug":
                        # update_Drug_xls(uploadFile, XlsSheet="Drug", upload=False, uploaduser=request.user, lower=True)
                    context["excel_upload_info"]="Saved in Drug Datatable"
                        
                elif process_name=="heatmap":
                    table= convert_heatmap(uploadFile, XlsSheet="Heatmap", upload=False, uploaduser=request.user, lower=True)
                    
                    context['table']=table
                    context["excel_upload_info"]="Convert data to heatmap"
                            
        except Exception as err:
            context["excel_upload_info"]=f"{err} cause failed import"
            # self.delete_file(file_name=filename)
        return render(request, self.template_name, context)
                  
