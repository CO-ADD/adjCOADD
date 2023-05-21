'''
View for uploading Vitek PDFs
'''
import concurrent
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
# get Progress 
def get_upload_progress(request):
    SessionKey = request.session.session_key
    progress =cache.get(SessionKey) or {'processed': 0, 'file_name':"",'total': 0, 'uploadpdf_version':1}   
    return JsonResponse(progress)
# 
# Cancel Thread
def cancel_upload(request):
    SessionKey = request.session.session_key
    # print(f"this is cancel process set{SessionKey}")
    cancel_flag_key = f'cancel_flag_{SessionKey}'
    cache.set(cancel_flag_key, True)
    cache.delete(request.session.session_key)
    # Reset/clear Wizard Session per request 
    del request.session['wizard_import__vitek_view']
    return HttpResponse('Upload canceled.')
# 
# Process View
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
        self.valLog=None
        if cache.get('vitek_process_step')==None:
            cache.set('vitek_process_step', 'upload_file')
            # self.storage.current_step =cache.get('vitek_process_step')
        
    def process_step(self, form):
        cache.set('vitek_process_step', self.steps.next)
        step=cache.get('vitek_process_step')
        print(f"next step: {step}")

        current_step = self.steps.current
        request = self.request
        # print(request.session.keys()) 
        # config session save and process cancel 
        SessionKey=self.request.session.session_key
        cancel_flag_key = f'cancel_flag_{SessionKey}'
        cache.set(cancel_flag_key, False)
        print(current_step)
        if current_step == 'upload_file':
            cache.set(f'current_step_{self.request.session.session_key}', current_step)
            DirName = file_location(instance=request.user)  # define file store path during file process
            files = []
            if form.is_valid():
                if 'upload_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('upload_file-multi_files'))          
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=DirName)
                    filename = fs.save(f.name, f)
                    self.filelist.append(filename)

                # Parse PDF and Validation
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(upload_VitekPDF_List, request, SessionKey, DirName, self.filelist, OrgBatchID=self.orgbatch_id, upload=False, appuser=request.user)
                    # block until the function finishes
                    self.valLog = future.result()  
                if self.valLog.nLogs['Error'] >0 :
                    dfLog = pd.DataFrame(self.valLog.get_aslist(logTypes= ['Error']))#convert result in a table
                    self.storage.extra_data['Confirm_to_Save'] = False
                else:
                    dfLog = pd.DataFrame(self.valLog.get_aslist())
                    self.storage.extra_data['Confirm_to_Save'] = True
                # self.storage.extra_data['valLog']=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False) 
                html=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False)
                cache.set(f'vitek_valLog_{request.user}', html)      
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['DirName'] = DirName          
            else:
                return render(request, 'ddrug/importhandler_vitek.html', context)

        elif current_step == 'step1': # recheck and save to DB
            # cache.set(f'current_step_{self.request.session.session_key}', current_step)
            # print("step validation again")
            upload=self.storage.extra_data['Confirm_to_Save']
            form = VitekValidation_StepForm(request.POST)
            self.organism_batch=request.POST.get("upload_file-orgbatch_id") #get organism_batch
            DirName=self.storage.extra_data['DirName'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name   

            with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(upload_VitekPDF_List, request, SessionKey, DirName, self.filelist, OrgBatchID=self.orgbatch_id, upload=upload, appuser=request.user)
                    # block until the function finishes
                    self.valLog = future.result()
           
            if self.valLog.nLogs['Error'] >0 :
                dfLog = pd.DataFrame(self.valLog.get_aslist(logTypes= ['Error']))#convert result in a table
            else:
                dfLog = pd.DataFrame(self.valLog.get_aslist())
            # self.storage.extra_data['valLog']=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False)
            html=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False)
            cache.set(f'vitek_valLog_{request.user}', html)  
                            
        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        print("Finalize")
        # Redirect to the desired page after finishing
        # delete uploaded files and clear cache
        filelist=self.storage.extra_data['filelist']
        for f in filelist:
            self.delete_file(f)
            print(filelist)
        cache.delete(self.request.session.session_key)
        cache.delete('vitek_process_step')
        cache.delete(f'vitek_valLog_{request.user}')
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        # save information to context,
        # then display in templates  
        context = super().get_context_data(form=form, **kwargs)
        context['step1']=self.name_step1
        current_step = self.steps.current     

        context['validation_result'] = cache.get(f'vitek_valLog_{self.request.user}') or "*/*"
        if current_step == 'step1':
            context['Confirm_to_Save']=self.storage.extra_data['Confirm_to_Save']
        # if current_step == 'finalize':
        #     context['validation_result'] = self.storage.extra_data['valLog']
        return context
    
    # def get(self, request, *args, **kwargs):
    #     # Check if a step was saved in the cache
    #     current_step =cache.get('vitek_process_step')
    #     print(f"current_step {current_step}")
    #     if current_step is not None:
    #         self.storage.current_step = current_step
    #     return super().get(request, *args, **kwargs)