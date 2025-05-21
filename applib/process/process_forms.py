import os
from django import forms
from django.shortcuts import HttpResponse, render, redirect
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict
from apputil.utils.views_base import SuperUserRequiredMixin, WriteUserRequiredMixin
from apputil.utils.files_upload import validate_file,file_location, OverwriteStorage


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
# 
#
    multi_files = SingleFileField(label='Select one file', 
                                  validators=[validate_file], 
                                  required=False)

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
                        validator(file)
                    except ValidationError as e:
                        self.add_error(field, f"{file.name}: {str(e)}")
                        
        if len(uploadfiles)<1: 
            self.add_error('single_file', "Select one file")
            raise forms.ValidationError("No files selected")
        return cleaned_data

# --------------------------------------------------------------------------------------------------
class SelectMultipleFiles_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = MultipleFileField(label='Select one or multiple files', 
                                  validators=[validate_file], 
                                  required=False)

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
            self.add_error('multi_files', "Select at least one file")
            raise forms.ValidationError("No files selected")
        return cleaned_data

# --------------------------------------------------------------------------------------------------
class Upload_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    confirm = forms.BooleanField(required=True, help_text="Confirm to upload Data")
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['confirm'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}

# --------------------------------------------------------------------------------------------------
class Finalize_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    pass

# =================================================================
# Process View
# -----------------------------------------------------------------
class Process_View(WriteUserRequiredMixin,SessionWizardView):
    
    COLUMN_FIELDS = ['Type','Note','Item','Filename','Help']
    
    process_step_names = ["Upload"]
    
    process_step1="Upload" # step label in template
    # define more steps name
    #... 
    # define each step's form
    
    process_form_list = [
        ('select_files', None),
        ('upload_files',None),
        # add more step -> StepForm
        ('finalize', None),
    ]
    # define template
    template_name = None

    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_list=[]
        self.file_dir=None
        self.valLog=None
        self.upload=False
        #self.html_columns = ['Type','Note','Item','Filename','Help']

    # File Processing and Validation and Upload 
    def file_process_handler(self, request, *args, **kwargs):
        #
        # request : HTML request
        #  file_dir : Folder for file(s)
        #  file_list: List of file names
        #  form_data: 
        #  upload : False - only validates the data in the files
        #           True - validates and uploads the data in the files
        #  appuser : logged in user 
        
        pass

    # ----------------------------------------------------
    # Main Process Step method
    def process_step(self, form):
    # ----------------------------------------------------
        current_step = self.steps.current
        request = self.request

        # Step 1 - Select File and Validate Data
        #---------------------------------------
        if current_step == 'select_files':
            context={}
            self.storage.extra_data['validation_result']="-"
            
            self.file_dir = file_location(instance=request.user)  # define file store path during file process
            _files = []

            if form.is_valid():
                if 'select_file-multi_files' in request.FILES:
                    _files.extend(request.FILES.getlist('select_file-multi_files'))        
                      
                # Get clean FileList
                for f in _files:
                    fs = OverwriteStorage(location=self.file_dir)
                    _filename = fs.save(f.name, f)
                    self.file_list.append(_filename)

                # Parse and Validate Files
                self.valLog=self.file_process_handler(request, 
                                                      self.file_dir, self.file_list, 
                                                      form_data=form.cleaned_data, 
                                                      upload=self.upload, appuser=request.user) 
                
                # Validation output Error
                if self.valLog.nLogs['Error'] >0 :
                    # Converts valLog Result into a table - no upload allowed
                    dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)
                    self.storage.extra_data['confirm_to_upload'] = False
                    
                # Validation output Warnings and Info - upload possible
                elif self.valLog.nLogs['Error'] <=0:
                    print(f"error is : {self.valLog.nLogs}")
                    try:
                        dfLog = self.valLog.get_ashtml(columns=self.html_columns)
                        self.storage.extra_data['confirm_to_upload'] = True
                    except Exception as err:
                        dfLog=f"{err}"
                        print(dfLog)
                        self.storage.extra_data['confirm_to_upload'] = False
                else:
                    dfLog = self.valLog.nLogs.get('Error') or 'No object exists, Is this a correct data file?'

                # Store Validation and File_Dir/List
                self.storage.extra_data['validation_result'] = dfLog
                self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) checked for errors." 
                self.storage.extra_data['file_list'] = self.file_list
                self.storage.extra_data['file_dir'] = self.file_dir          
            else:
                self.storage.extra_data['validation_result']="No files selected"
                return render(request, self.template_name, context)

        # Step 2 - Upload Data 
        #---------------------------------------
        elif current_step == 'upload_files': # recheck and save to DB
            form =self.form_list['upload_files'](request.POST)
            
            self.upload=True
            self.file_dir=self.storage.extra_data['file_dir'] #get file path
            self.file_list=self.storage.extra_data['file_list'] #get files' name 
            
            # Parse, Validate and Upload Files
            self.valLog=self.file_process_handler(request, 
                                                  self.file_dir, self.file_list, 
                                                  form_data=request.POST, 
                                                  upload=self.upload, appuser=request.user)
            
            if self.valLog.nLogs['Error'] >0 :
                dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)#convert result in a table
            else:
                dfLog = self.valLog.get_ashtml(columns=self.html_columns)

            # Store Validation and File_Dir/List
            self.storage.extra_data['validation_result'] = dfLog  
            self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) Uploaded." 

        return self.get_form_step_data(form)