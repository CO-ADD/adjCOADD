import os
import shutil
from django import forms
from django.shortcuts import HttpResponse, render, redirect
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect

from applib.process.process_stepforms import (MultipleFileInput, SingleFileInput,
                                              MultipleFileField, SingleFileField,
                                              SelectSingleFile_StepForm, SelectMultipleFiles_StepForm, SelectSingleFileFolder_StepForm,
                                              Upload_StepForm, Generate_StepForm,
                                              Finalize_StepForm)
  
from applib.django.base.views import SuperUserRequiredMixin, WriteUserRequiredMixin
from apputil.utils.files_upload import validate_file, file_location, OverwriteStorage
#from apputil.utils.validation_log import Validation_Log
from applib.logging.validation_log import Validation_Log




# =================================================================
# Progress Utilities
# -----------------------------------------------------------------
def Start_Progress(request):
    return render(request, 'modal/progress_partial_modal.html')

def Update_Progress(request,):
    current_progress = 50 # Replace with actual logic
    total_progress = 100
    return render(request, 'modal/progress_partial_modal.html', {'current': current_progress, 'total':total_progress})



# --------------------------------------------------------------------------------------------------
class Process_View(WriteUserRequiredMixin,SessionWizardView):
# --------------------------------------------------------------------------------------------------

    process_name = 'File Upload'
    
    name_step1="Upload" # step label in template
    # define more steps name
    #... 
    # define each step's form
    form_list = [
        ('select_file', None),
        ('upload',None),
        # add more step -> StepForm
        ('finalize', None),
    ]
    # define template
    template_name = None
    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    # ----------------------------------------------------
    def __init__(self, *args, **kwargs):
        print(f" [Process_View.__init__] ")
        super().__init__(*args, **kwargs)
        self.file_list={}
        self.file_dir=None
        self.pk = None
        self.valLog=None
        self.upload=False
        self.overwrite=False
        self.html_columns = Validation_Log.LOG_FIELDS
    
    # ----------------------------------------------------
    def get_object(self):
        self.pk = self.kwargs.get('pk')

        print(f" [get_object] {self.pk}")
        self.object = get_object_or_404(self.model, pk=self.pk)

    # ----------------------------------------------------
    def file_process_handler(self, request, *args, **kwargs):
        print(f" [file_process_handler] not implemented")
        
    # ----------------------------------------------------
    def file_process_finalizer(self, request, *args, **kwargs):
        print(f" [file_process_finalizer] not implemented")

    # ----------------------------------------------------
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request

        print(f" [process_step] Step: {current_step} Request: {request} ")
        
        # First Step - Select File(s) -> self.filelist[{file_field as per select form}]
        if current_step == 'select_file':
            context={}
            self.storage.extra_data['validation_result']="-"
            self.file_dir = file_location(instance=request.user)  # define file store path during file process

            _files = {}
            
            if 'object_pk' in self.storage.extra_data:
                # ProcessView for Existing Entry
                self.pk = self.storage.extra_data['object_pk']
            else:
                # ProcessView for New Entry
                self.pk = None
                
            if form.is_valid():
                print(f" [Process_View.process_step] Valid Form {self.pk}")
                
                # Get list of files for each select_file-{file-field}
                for _key in request.FILES:
                    _field = None
                    if 'select_file' in _key:
                        _field = _key.replace('select_file-','')
                    
                    if _field:
                        _files[_field] = request.FILES.getlist(_key)
                        
                #print(f" [process_step] _files {_files}")
                # if 'select_file-multi_files' in request.FILES:
                #     files.extend(request.FILES.getlist('select_file-multi_files'))
                              
                # Get clean/save files and create file_list[{file-field}] 
                self.file_list = {}
                for _field in _files:
                    self.file_list[_field] = []
                    for f in _files[_field]:
                        fs = OverwriteStorage(location=self.file_dir)
                        filename = fs.save(f.name, f)
                        self.file_list[_field].append(filename)

                print(f" [process_step] file_list {self.file_list}")
                # Parse and Validation
                self.valLog=self.file_process_handler(request, 
                                                      self.file_dir, self.file_list, 
                                                      form_data=form.cleaned_data, 
                                                      upload=self.upload, appuser=request.user) 
                
                #convert result in a table
                if self.valLog.n_logs['Error'] >0 :
                    #dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)
                    dfLog = self.valLog.get_ashtml(columns=self.html_columns)
                    self.storage.extra_data['confirm_to_upload'] = False
                    
                elif self.valLog.n_logs['Error'] <=0:
                    #print(f" [{self.process_name}] valLog: {self.valLog.n_logs}")
                    # try:
                    dfLog = self.valLog.get_ashtml(columns=self.html_columns)
                    self.storage.extra_data['confirm_to_upload'] = True
                    # except Exception as err:
                    #     dfLog=f"{err}"
                    #     self.storage.extra_data['confirm_to_upload'] = False
                else:
                    dfLog = self.valLog.nLogs.get('Error') or 'No object exists, Is this a correct data file?'

                self.storage.extra_data['validation_result'] = dfLog
                self.storage.extra_data['validation_message']= f" File(s) checked for errors: {len(self.file_list)}" 
                #self.storage.extra_data['validation_help']= self.message_html['upload']
                self.storage.extra_data['file_list'] = self.file_list
                self.storage.extra_data['file_dir'] = self.file_dir
            else:
                self.storage.extra_data['validation_result']="No files selected"
                #print(" [Process_View.process_step] Not Valid")
                return render(request, self.template_name, context)

        elif current_step == 'upload': # recheck and save to DB
            #print(f" [Process_View.process_step] Request.Post: {request.POST}")
            #print(f" [Process_View.process_step] Request.Post: {request.POST}")
            if form.is_valid():
                #print(f" [Process_View.process_step] Form CleanedData: {form.cleaned_data}")

                self.file_dir=self.storage.extra_data['file_dir'] #get file path
                self.file_list=self.storage.extra_data['file_list'] #get files' name  
                self.pk = self.storage.extra_data['object_pk']
                
                #print(f" [Process_View.process_step] {current_step} PK: {self.pk}")
                #print(f" [Process_View.process_step] {current_step} Form: {form}")

                self.valLog=self.file_process_handler(request, 
                                                    self.file_dir, self.file_list, 
                                                    form_data=form.cleaned_data, 
                                                    upload=self.upload, appuser=request.user)
                
                if self.valLog.n_logs['Error'] >0 :
                    dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)
                else:
                    dfLog = self.valLog.get_ashtml(columns=self.html_columns)

                self.storage.extra_data['validation_result'] = dfLog  
                self.storage.extra_data['validation_message']= f" File(s) uploaded: {len(self.file_list)} "
                #self.storage.extra_data['validation_help'] = "Done" 

        elif current_step == 'finalize':
            # self.file_dir=self.storage.extra_data['file_dir'] #get file path
            # self.file_list=self.storage.extra_data['file_list'] #get files' name  
            self.pk = self.storage.extra_data['object_pk']
            print(f" [process_step] {current_step} PK: {self.pk}")
            self.file_process_finalizer(request, self.pk)
                            
        return self.get_form_step_data(form)

    # ----------------------------------------------------
    def done(self, form_list, **kwargs):
        #print(f" [Process_View.done] ")
        # Redirect to the desired page after finishing
        file_dir=self.storage.extra_data['file_dir']
        #print(file_dir)
        if file_dir:
            try:
                shutil.rmtree(file_dir)
                
            except FileNotFoundError as err:
                print(err)
            except Exception as err:
                print(err)
        return redirect(self.request.META['HTTP_REFERER'])


    # ----------------------------------------------------
    def get_context_data(self, form, **kwargs):
        context = super().get_context_data(form=form, **kwargs)
        
        self.pk = self.kwargs.get('pk',None)
        if self.pk:
            self.storage.extra_data['object_pk'] = self.pk
        
        #print(f" [Process_View.get_context_data] PK: {self.pk} ")
        # save information to context,
        # then display in templates
          
        #context['step1']=self.name_step1
        current_step = self.steps.current
        context['validation_message'] = self.storage.extra_data.get('validation_message', None)
        context['help_message'] = self.storage.extra_data.get('help_message', None)

        #print(f" [Process_View.get_context_data] Step: [{current_step}] for PK:{self.pk} ")
        
        context['pk'] = self.storage.extra_data.get('object_pk', None)
        if current_step == 'select_file':
            context['validation_result']=""
            #context['pk'] = self.pk
        elif current_step == 'upload_file':
            context['validation_result']=""
            #context['pk'] = self.storage.extra_data.get('object_pk', None)
        else:
            context['validation_result'] = self.storage.extra_data.get('validation_result', None)
            context['confirm_to_upload'] = self.storage.extra_data.get('confirm_to_upload', None)
            context['help_text']=""
            #context['pk'] = self.storage.extra_data.get('object_pk', None)
            
        #print(f" [Process_View.get_context_data] Step: [{current_step}] Validation_result: {context['validation_result']}")
        return context
    
    

