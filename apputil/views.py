import os
from formtools.wizard.views import SessionWizardView
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
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from dorganism.models import Organism, Taxonomy
from ddrug.models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
from dgene.models import Gene, WGS_CheckM, WGS_FastQC, ID_Pub, ID_Sequence

from apputil.forms import AppUserfilter, Dictionaryfilter, ApplicationUser_form, Dictionary_form, Login_form, Image_form
from apputil.models import ApplicationUser, Dictionary, Image, Document
from apputil.utils.views_base import SuperUserRequiredMixin, permission_not_granted, SimplecreateView, HtmxupdateView
from apputil.utils.filters_base import FilteredListView
from apputil.utils.files_upload import Importhandler

## =================================APP Home========================================

from django.http import JsonResponse
from django.core.files.storage import default_storage

@login_required
def createImage(request, pk):
    object_ = get_object_or_404(Organism, organism_id=pk)
    kwargs={}
    kwargs['user']=request.user
    form=Image_form 

    # Handle AJAX file upload
    if request.method == 'POST' and request.headers.get('x-requested-with') == 'XMLHttpRequest':
        image_file = request.FILES.get('image_file')
        if image_file:
            file_path = default_storage.save(image_file.name, image_file)
            file_name = os.path.basename(file_path)
            file_type = image_file.content_type
            response = {
                'name': file_name,
                'type': file_type,
                'path': file_path,
            }
            return JsonResponse(response)

    # Handle form submission
    elif request.method == 'POST':
        form=Image_form(request.POST, request.FILES)                
        if form.is_valid():    
            instance=form.save(commit=False)
            if instance.image_file:
                instance.save(**kwargs)
                print(f'saved:{instance.image_file}')
                object_.assoc_images.add(instance)
                object_.save(**kwargs)
            else:
                messages.warning(request, 'No image file.')
            return redirect(request.META['HTTP_REFERER'])
        else:
            messages.warning(request, f'Update failed due to {form.errors} error')
            return redirect(request.META['HTTP_REFERER'])

    context={
        "image_form":form,
        "object":object_,
    }

    return render(request, "apputil/addimg.html", context)
# from django import forms
# from django.core.files.storage import DefaultStorage
# class Image_uploadform(forms.ModelForm):
#     image_file = forms.ImageField(label='Select an image', 
#                                 #   validators=[validate_file], 
#                                   required=True)
#     class Meta:
#         model=Image
#         fields=['image_file']

# class Image_infoform(forms.ModelForm):
   
#     class Meta:
#         model=Image
#         exclude=['image_file']


# class SimpleUpload_WizardView(SessionWizardView):
#     form_list=[Image_uploadform, Image_infoform]
#     template_name='utils/simpleupload.html'
#     file_storage = DefaultStorage()
    
#     def get(self, request, *args, **kwargs):
#         self.object = get_object_or_404(Organism, organism_id=self.kwargs['pk'])  # get your model instance here
#         return super().get(request, *args, **kwargs)

#     def get_context_data(self, form, **kwargs):
#         context = super().get_context_data(form=form, **kwargs)
#         context.update({'object': self.object})  # add your model instance to the context
#         return context

#     def get_form_initial(self, step):
#         initial = super().get_form_initial(step)
#         if step == '1':  # the step index of Image_infoform
#             image_data = self.get_cleaned_data_for_step('0')  # get image data from Image_uploadform
#             if image_data:
#                 initial.update({
#                     'image_type': image_data['image_file'].content_type, 
#                     'image_name': image_data['image_file'].name,
#                     'image_source': image_data['image_file'].name,
#                     'image_desc': image_data['image_file'].name,
                  
#                 })
#         return initial

#     def done(self, form_list, **kwargs):
#         image_upload_form_data = self.get_cleaned_data_for_step('0')  # get data from Image_uploadform
#         image_info_form = form_list[-1]  # this should be Image_infoform
#         if image_info_form.is_valid():
#             image_info_data = image_info_form.cleaned_data
#             image = Image()  # create new Image instance
#             image.image_file = image_upload_form_data['image_file']
#         for field, value in image_info_data.items():
#             setattr(image, field, value)  # set other fields from Image_infoform data
#         image.save()  # save the image
#         related_model=get_object_or_404(Organism, organism_id=self.kwargs['pk'])
#         related_model.assoc_images.add(image)

#         return HttpResponse("form submitted!")

@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def deleteImage(req, pk):
    print('deleting...')
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(Image, id=pk)
    try:
        object_.delete(**kwargs)
        print('deleted')
    except Exception as err:
        print(err)
    return redirect(req.META['HTTP_REFERER'])

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

# =========================== Export CSV View =============================
from .utils.views_base import DataExportBaseView
class DataExportView(DataExportBaseView):
    pass

# =============Import Dictionary and appUsers via Excel==================
from .utils.files_upload import FileUploadForm, OverwriteStorage, file_location
from .utils.validation_log import Validation_Log
from .utils.data_visual import convert_heatmap
from .utils.form_wizard_tools import SelectFile_StepForm

class Importhandler_apputils(Importhandler):
    form_class=SelectFile_StepForm
    template_name='apputil/importhandler_excel.html'
    upload_model_type=None
    lObject={}   # store results parsed by all uploaded excel files with key-filename, value-parsed result array
    # vLog = Validation_Log("Vitek-pdf")
    
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
                  
