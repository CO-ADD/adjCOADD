import os
import shutil
from django import forms
from django.shortcuts import HttpResponse, render, redirect
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect

from applib.django.base.views import SuperUserRequiredMixin, WriteUserRequiredMixin
from apputil.utils.files_upload import validate_file, file_location, OverwriteStorage
#from apputil.utils.validation_log import Validation_Log
from applib.logging.validation_log import Validation_Log


# =================================================================
# Utilities Forms
# -----------------------------------------------------------------
class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True

# -----------------------------------------------------------------
class SingleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = False

# -----------------------------------------------------------------
class MultipleFileField(forms.FileField):
    
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = single_file_clean(data, initial)
        return result

# -----------------------------------------------------------------
class SingleFileField(forms.FileField):
    
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("widget", SingleFileInput())
        super().__init__(*args, **kwargs)
       
        
    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            result = [single_file_clean(d, initial) for d in data]
        else:
            result = single_file_clean(data, initial)
        return result
    
# =================================================================
# Process Step Forms
# -----------------------------------------------------------------
class SelectSingleFile_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = SingleFileField(label='Upload File', 
                                  validators=[validate_file], 
                                  required=False,
                                  help_text="Select a single file")
        
    def clean(self):
        cleaned_data = super().clean()
        uploadfiles=[]
        # List of file fields to validate
        file_fields = ['multi_files',]
        
        # check if filelist is MultiValueDict
        if not isinstance(self.files, MultiValueDict):
            return cleaned_data
        
        for field in file_fields:
            files = self.files.getlist(f'select_file-{field}')
            
            uploadfiles.extend(files)

            for file in files:
                for validator in self.fields[field].validators:
                    try:
                        validator(file)
                    except ValidationError as e:
                        self.add_error(field, f"{file.name}: {str(e)}")
                        
        if len(uploadfiles)<1: 
            self.add_error('multi_files', "No file selected. Please select a file")
            raise forms.ValidationError("No files selected")
        
        return cleaned_data

# --------------------------------------------------------------------------------------------------
class SelectMultipleFiles_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = MultipleFileField(label='Upload Files', 
                                  validators=[validate_file], 
                                  required=False,
                                  help_text="Select one or multiple files")

    def clean(self):
        cleaned_data = super().clean()
        uploadfiles=[]
        # List of file fields to validate
        file_fields = ['multi_files',]
        
        # check if filelist is MultiValueDict
        if not isinstance(self.files, MultiValueDict):
            return cleaned_data
        
        for field in file_fields:
            files = self.files.getlist(f'select_file-{field}')
            
            uploadfiles.extend(files)

            for file in files:
                for validator in self.fields[field].validators:
                    try:
                        validator(file)
                    except ValidationError as e:
                        self.add_error(field, f"{file.name}: {str(e)}")
                        
        if len(uploadfiles)<1: 
            self.add_error('multi_files', "No file(s) selected. Please select at least one file")
            raise forms.ValidationError("No files selected")
        return cleaned_data

# --------------------------------------------------------------------------------------------------
class Upload_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    upload = forms.BooleanField(initial=False, required=False, help_text="Upload New Data to Database")
    overwrite = forms.BooleanField(initial=False, required=False, help_text="Overwrite Existing Data as well")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['upload'].label = "Upload New Data"
        self.fields['overwrite'].label = "Overwrite Existing Data"
        # self.fields['upload'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}
        # self.fields['overwrite'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}

# --------------------------------------------------------------------------------------------------
class Finalize_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    pass


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
        #print(f" [Process_View.__init__] ")
        super().__init__(*args, **kwargs)
        self.file_list=[]
        self.file_dir=None
        self.pk = None
        self.valLog=None
        self.upload=False
        self.overwrite=False
        self.html_columns = Validation_Log.LOG_FIELDS
    
    # ----------------------------------------------------
    def get_object(self):
        #print(f" [Process_View.get_object] ")
        self.pk = self.kwargs.get('pk')
        self.object = get_object_or_404(self.model, pk=self.pk)

    # ----------------------------------------------------
    def file_process_handler(self, request, *args, **kwargs):
        pass
        
    # ----------------------------------------------------
    def file_process_finalizer(self, request, *args, **kwargs):
        pass

    # ----------------------------------------------------
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request

        #print(f" [Process_View.process_step] Step: {current_step} Request: {request} ")

        if current_step == 'select_file':
            context={}
            self.storage.extra_data['validation_result']="-"
            self.file_dir = file_location(instance=request.user)  # define file store path during file process
            files = []
            
            self.pk = self.storage.extra_data['object_pk']
            if form.is_valid():
                #print(f" [Process_View.process_step] Valid Form {self.pk}")

                if 'select_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('select_file-multi_files'))
                              
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=self.file_dir)
                    filename = fs.save(f.name, f)
                    self.file_list.append(filename)

                # Parse and Validation
                self.valLog=self.file_process_handler(request, 
                                                      self.file_dir, self.file_list, 
                                                      form_data=form.cleaned_data, 
                                                      upload=self.upload, appuser=request.user) 
                
                #convert result in a table
                #if self.valLog.nLogs['Error'] >0 :
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
            #print(f" [Process_View.process_step] {current_step} PK: {self.pk}")
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
    
    
# # =================================================================
# # Process View
# # -----------------------------------------------------------------
# class XProcess_View(WriteUserRequiredMixin,SessionWizardView):
    
#     COLUMN_FIELDS = ['Type','Note','Item','Filename','Help']
    
#     STEPS = {
#         '01_Select_Files' : {'form':None,'name': 'Select files'},
#         '02_Upload_Files' : {'form':None,'name': 'Upload files'},
#         '03_Finalize'     : {'form':None,'name': 'Finish'}
#     }

#     file_storage = FileSystemStorage(location='/tmp/')
    
#     # -----------------------------------------------
#     model = None
#     template_name = None

#     process_step_names = ["Upload"]
    
#     process_step1="Upload" # step label in template
#     # define more steps name
#     #... 
#     # define each step's form
    
#     process_form_list = [
#         ('select_files', None),
#         ('upload_files',None),
#         # add more step -> StepForm
#         ('finalize', None),
#         ]
    

#     # Define a file storage for handling file uploads
#     # ----------------------------------------------------
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.file_list=[]
#         self.file_dir=None
#         self.valLog=None
#         self.upload=False
#         #self.html_columns = ['Type','Note','Item','Filename','Help']
        
#     # ----------------------------------------------------
#     def get_object(self, queryset=None):
#         self.pk = self.kwargs.get('pk')
#         self.object = get_object_or_404(self.model, pk=self.pk)

#     # ----------------------------------------------------
#     # File Processing and Validation and Upload 
#     def file_process_handler(self, request, *args, **kwargs):
#         #
#         # request : HTML request
#         #  file_dir : Folder for file(s)
#         #  file_list: List of file names
#         #  form_data: 
#         #  upload : False - only validates the data in the files
#         #           True - validates and uploads the data in the files
#         #  appuser : logged in user 
        
#         pass

#     # ----------------------------------------------------
#     # Main Process Step method
#     def process_step(self, form):
#     # ----------------------------------------------------
#         current_step = self.steps.current
#         request = self.request

#         print(f'[Process_View] {current_step} ')
        
#         # Step 1 - Select File and Validate Data
#         #---------------------------------------
#         if current_step == '01_Select_Files':
            
#             context={}
#             self.storage.extra_data['validation_result']="-"
            
#             self.file_dir = file_location(instance=request.user)  # define file store path during file process
#             _files = []

#             if form.is_valid():
#                 if 'select_file-multi_files' in request.FILES:
#                     _files.extend(request.FILES.getlist('select_file-multi_files'))        
                      
#                 # Get clean FileList
#                 for f in _files:
#                     fs = OverwriteStorage(location=self.file_dir)
#                     _filename = fs.save(f.name, f)
#                     self.file_list.append(_filename)

#                 # Parse and Validate Files
#                 self.valLog=self.file_process_handler(request, 
#                                                       self.file_dir, self.file_list, 
#                                                       form_data=form.cleaned_data, 
#                                                       upload=self.upload, appuser=request.user) 
#                 print(f" [Process_Step] valLog : {self.valLog.nLogs}") 
#                 # Validation output Error
#                 if self.valLog.nLogs['Error'] >0 :
#                     # Converts valLog Result into a table - no upload allowed
#                     dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)
#                     self.storage.extra_data['confirm_to_upload'] = False
                    
#                 # Validation output Warnings and Info - upload possible
#                 elif self.valLog.nLogs['Error'] <=0:
#                     try:
#                         dfLog = self.valLog.get_ashtml(columns=self.html_columns)
#                         self.storage.extra_data['confirm_to_upload'] = True
#                     except Exception as err:
#                         dfLog=f"{err}"
#                         print(dfLog)
#                         self.storage.extra_data['confirm_to_upload'] = False
#                 else:
#                     dfLog = self.valLog.nLogs.get('Error') or 'No object exists, Is this a correct data file?'

#                 # Store Validation and File_Dir/List
#                 self.storage.extra_data['validation_result'] = dfLog
#                 self.storage.extra_data['validation_message']= f" {len(self.file_list)} file(s) checked for errors." 
#                 self.storage.extra_data['file_list'] = self.file_list
#                 self.storage.extra_data['file_dir'] = self.file_dir          
#             else:
#                 self.storage.extra_data['validation_result']="No files selected"
#                 return render(request, self.template_name, context)

#         # Step 2 - Upload Data 
#         #---------------------------------------
#         elif current_step == '02_Upload_Files': # recheck and save to DB
#             form =self.form_list['upload_files'](request.POST)
            
#             self.upload=True
#             self.file_dir=self.storage.extra_data['file_dir'] #get file path
#             self.file_list=self.storage.extra_data['file_list'] #get files' name 
            
#             # Parse, Validate and Upload Files
#             self.valLog=self.file_process_handler(request, 
#                                                   self.file_dir, self.file_list, 
#                                                   form_data=request.POST, 
#                                                   upload=self.upload, appuser=request.user)
            
#             if self.valLog.nLogs['Error'] >0 :
#                 dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)#convert result in a table
#             else:
#                 dfLog = self.valLog.get_ashtml(columns=self.html_columns)

#             # Store Validation and File_Dir/List
#             self.storage.extra_data['validation_result'] = dfLog  
#             self.storage.extra_data['validation_message']= f" {len(self.file_list)} file(s) Uploaded." 

#         return self.get_form_step_data(form)

#     # ----------------------------------------------------
#     def done(self, form_list, **kwargs):
#     # ----------------------------------------------------
#         # Redirect to the desired page after finishing
#         file_dir=self.storage.extra_data['file_dir']
#         if file_dir:
#             try:
#                 shutil.rmtree(file_dir)
                
#             except FileNotFoundError as err:
#                 print(err)
#             except Exception as err:
#                 print(err)
#         return redirect(self.request.META['HTTP_REFERER'])

    
#     # ----------------------------------------------------
#     def get_context_data(self, form, **kwargs):
#     # ----------------------------------------------------
#         context = super().get_context_data(form=form, **kwargs)
#         # save information to context,
#         # then display in templates  
#         context['step1']=self.name_step1
#         current_step = self.steps.current
#         context['validation_message'] = self.storage.extra_data.get('validation_message', None)
#         if current_step == '01_Upload_File':
#             context['validation_result']="Select Files"
#         else:
#             context['validation_result'] = self.storage.extra_data.get('validation_result', None)
#             context['confirm_to_upload']=self.storage.extra_data.get('confirm_to_upload', None)
#         print(f"[ImportHandler_View] {current_step} validation_result: {context['validation_result']}")
#         return context
