'''
View for uploading Vitek PDFs
'''
from asgiref.sync import async_to_sync, sync_to_async
import pandas as pd
import threading

from django import forms
from django.conf import settings
from django.core.cache import cache
from django.core.files.storage import FileSystemStorage
from django.http import JsonResponse
from django.shortcuts import HttpResponse, render, redirect

from apputil.utils.files_upload import file_location, OverwriteStorage
from apputil.utils.form_wizard_tools import ImportHandler_WizardView, UploadFileForm, StepForm_1, FinalizeForm
from apputil.utils.validation_log import * 
from dorganism.models import Organism_Batch
from ddrug.models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub
from ddrug.utils.vitek import upload_VitekPDF_List

#==  VITEK Import View =============================================================

# customized Form
class VitekValidation_StepForm(StepForm_1):
    orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',
                                       widget=forms.Select(attrs={'class': 'form-control'}), 
                                       required=False, 
                                       help_text='(**Optional)',
                                       queryset=Organism_Batch.objects.filter(astatus__gte=0), )
    field_order = ['orgbatch_id', 'confirm']
# 

# 
def get_upload_progress(request):
    # session_key = f'upload_progress_{request.user}'
    SessionKey = request.session.session_key
    progress =cache.get(SessionKey) or {'processed': 0, 'file_name':"",'total': 0, 'uploadpdf_version':0}   
    return JsonResponse(progress)
# 
# For get result
# 
def fetchResult(request):
    cache_key = f'valLog_{request.user}'
    valLog = cache.get(cache_key) or None
    if valLog is not None:
        return JsonResponse({
            'results_ready': True,
            'validation_result': valLog,
        })
    else:
        return JsonResponse({'results_ready': False})
# 
# call uploading function(upload_Vitek_List ....)
def uplod_utils(Request, SessionKey, DirName,FileList,OrgBatchID=None,upload=False,appuser=None):
    # Cancel Process, for restart a uploading
    cancel_flag_key = f'cancel_flag_{SessionKey}' # Setting up the cancel flag key
    if cache.get(cancel_flag_key):
        cache.delete(f'valLog_{self.request.user}')
        cache.delete(SessionKey)
        return None
    valLog=upload_VitekPDF_List(Request, SessionKey, DirName,FileList,OrgBatchID=None,upload=False,appuser=None)
    print(valLog)
    # Create vLog save to Cache
    cache_key = f'valLog_{Request.user}'
    if valLog.nLogs['Error'] >0 :
        dfLog = pd.DataFrame(valLog.get_aslist(logTypes= ['Error']))#convert result in a table
        Confirm_to_Save = False
    else:
        dfLog = pd.DataFrame(valLog.get_aslist())
        Confirm_to_Save = True
    valLog=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "overflow-auto"], index=False)
    cache.set(cache_key, {'Confirm_to_Save':Confirm_to_Save, 'valLog':valLog})
    return None
# 
# Cancel Thread
def cancel_upload(request):
    SessionKey = request.session.session_key
    print(f"this is cancel process set{SessionKey}")
    cancel_flag_key = f'cancel_flag_{SessionKey}'
    cache.set(cancel_flag_key, True)
    cache.delete(f'valLog_{request.user}')
    cache.delete(request.session.session_key)
    # Reset/clear Wizard Session per request 
    del request.session['wizard_import__vitek_view']
    return HttpResponse('Upload canceled.')
# 
class Import_VitekView(ImportHandler_WizardView):
    
    name_step1="Validation" # step label in template
    # define more steps name
    #... 
    # define each step's form
    form_list = [
        ('upload_file', UploadFileForm),
        ('step1', VitekValidation_StepForm),
        # add more step -> StepForm
        ('finalize', FinalizeForm),
    ]
    # define template
    template_name = 'ddrug/importhandler_vitek.html'
    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filelist=[]
        self.orgbatch_id=None
        
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request
        print(request.session.keys()) 
        # config session save and process cancel 
        SessionKey=self.request.session.session_key
        cancel_flag_key = f'cancel_flag_{SessionKey}'
        cache.set(cancel_flag_key, False)
        
        if current_step == 'upload_file':
            # print(self.storage)
            print("upload and parse")
            # print(self.request.session.get('formtools_wizard'))
            DirName = file_location(instance=request.user)  # define file store path during file process
            files = []
            result_table=[]
            if form.is_valid():
                if 'upload_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('upload_file-multi_files'))          
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=DirName)
                    filename = fs.save(f.name, f)
                    self.filelist.append(filename)

                # Parse PDF 
                thread = threading.Thread(target=uplod_utils, args=(request, SessionKey, DirName, self.filelist ), kwargs={'OrgBatchID': self.orgbatch_id, 'upload': False, 'appuser':request.user})
                thread.start()
                print(f"Currently running threads: {threading.enumerate()}")             
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['DirName'] = DirName
                print(self.storage.extra_data['filelist'])          
            else:
                return render(request, 'ddrug/importhandler_vitek.html', context)

        elif current_step == 'step1': # first validation
            upload=False
            print("step validation again")
            form = VitekValidation_StepForm(request.POST)
            if form.is_valid():
                confirm = form.cleaned_data.get('confirm')
                # print(confirm)
                if confirm:
                    upload=confirm
                    print("data will be uploaded")
            self.organism_batch=request.POST.get("upload_file-orgbatch_id") #get organism_batch
            result_table=[] # first validation result tables
            DirName=self.storage.extra_data['DirName'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name   
            thread = threading.Thread(target=uplod_utils, args=(request, SessionKey, DirName, self.filelist ), kwargs={'OrgBatchID': self.orgbatch_id, 'upload':upload, 'appuser':request.user})
            thread.start()
            # print(f"Currently running threads: {threading.enumerate()}")                 
        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        print("Finalize")
        # Redirect to the desired page after finishing
        # delete uploaded files
        filelist=self.storage.extra_data['filelist']
        for f in filelist:
            self.delete_file(f)
            print(filelist)
        cache.delete(f'valLog_{self.request.user}')
        cache.delete(self.request.session.session_key)
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        # save information to context,
        # then display in templates  
        context = super().get_context_data(form=form, **kwargs)
        context['step1']=self.name_step1
        current_step = self.steps.current           
        return context
