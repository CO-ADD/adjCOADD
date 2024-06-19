
import os,shutil

from django import forms
from django.utils.datastructures import MultiValueDict
from django.core.exceptions import ValidationError
from django.core.files.storage import FileSystemStorage
from django.shortcuts import render,redirect
from formtools.wizard.views import SessionWizardView

from apputil.utils.validation_log import Validation_Log
from apputil.utils.views_base import SuperUserRequiredMixin, WriteUserRequiredMixin
from apputil.utils.files_upload import FileValidator,file_location, OverwriteStorage

"""
    Steps process driven by form submission requests. 
        For uploading data from Files (PDF, Excel)

        - 1: Select Files to upload -> parsing, validating
        - 2: Upload data -> parsing, validating and uploading to database
        - 3: Outcome -> show uploading outcome 

    Process variables in sessionView storages: 
        valLog               validation result table store in "validation-result"; 
        confirm-to-upload    enable upload data store in "validation-result";
        dirname              files store path store in "dirname";
        filelist             uploaded files' namelist store in "filelist";
        validation_message   summarize progress store in "validation_message";
        step_1               store step name.
"""

# ==================================================================================
# Step 1 - File upload - Form Fields
# ==================================================================================

#-----------------------------------------------------------------------------------
class MultipleFileInput(forms.ClearableFileInput):
#-----------------------------------------------------------------------------------
    allow_multiple_selected = True

#-----------------------------------------------------------------------------------
class SingleFileInput(forms.ClearableFileInput):
#-----------------------------------------------------------------------------------
    allow_multiple_selected = False

#-----------------------------------------------------------------------------------
class MultipleFileField(forms.FileField):
#-----------------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------------
class SingleFileField(forms.FileField):
#-----------------------------------------------------------------------------------
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

# ==================================================================================
# Step 1 - File upload Forms
# ==================================================================================
Data_FileValidator = FileValidator(#max_size=1024 * 100, 
                             content_types=('text/csv', 
                                            'application/pdf',
                                            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', 
                                            ))


# --------------------------------------------------------------------------------------------------
class SelectSingleFile_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = SingleFileField(label='Select one file', 
                                  validators=[Data_FileValidator], 
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
            self.add_error('single_file', "Select one file")
            raise forms.ValidationError("No files selected")
        return cleaned_data

# --------------------------------------------------------------------------------------------------
class SelectMultipleFiles_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = MultipleFileField(label='Select one or multiple files', 
                                  validators=[Data_FileValidator], 
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

# ==================================================================================
# Upload Handler
# ==================================================================================

# --------------------------------------------------------------------------------------------------
class UploadHandler_View(WriteUserRequiredMixin,SessionWizardView):
# --------------------------------------------------------------------------------------------------

    # Define self.variables - on startup, before any UploadHandler_View.__init__()
    file_storage = FileSystemStorage(location='/tmp/')
    template_name = None
    form_list = [
                        ('select_file', None),
                        ('upload',None),
                        # add more step -> StepForm
                        ('finalize', None),
                    ]
    # -----------------------------------------------------------------------------

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filelist=[]
        self.dirname=None
        self.valLog=None
        self.upload=False
        self.html_columns = ['Type','Note','Item','Filename','Help']
 
    # Customised for each upload type  
    def file_process_handler(self, request, *args, **kwargs) -> Validation_Log:
        pass
        
    # Processing Step 1 and 2     
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request

        # Step 1 - select_file
        if current_step == 'select_file':
            context={}
            self.storage.extra_data['validation_result']="-"

            self.dirname = file_location(instance=request.user)  # define file store path during file process
            _Files = []

            _confirm_to_upload = False
            _remove_files = False

            if form.is_valid():
                if 'select_file-multi_files' in request.FILES:
                    _Files.extend(request.FILES.getlist('select_file-multi_files'))

                # Get clean FileList
                for f in _Files:
                    fs = OverwriteStorage(location=self.dirname)
                    _Filename = fs.save(f.name, f)
                    self.filelist.append(_Filename)

                # -- Run file processing function
                self.valLog=self.file_process_handler(request, self.dirname, self.filelist, form_data=form.cleaned_data, upload=self.upload, appuser=request.user)

                # -- Convert valLog results from file processing into HTML table
                if self.valLog.nLogs['Error'] <=0:
                    _htmlLog = self.valLog.get_ashtml(columns=self.html_columns)
                    _confirm_to_upload = True
                elif self.valLog.nLogs['Error'] >0 :
                    _htmlLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)
                    _remove_files = True 
                else:
                    _htmlLog = self.valLog.nLogs.get('Error') or 'No data found. Is this a correct file ?'
                    _remove_files = True 

                self.storage.extra_data['validation_result'] = _htmlLog
                self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) checked for upload." 
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['dirname'] = self.dirname
                self.storage.extra_data['confirm_to_upload'] = _confirm_to_upload

                # Clean up files in case of Errors
                if _remove_files:
                    self._remove_uploaded_files()          
          
            else:
                self.storage.extra_data['validation_result']= "No files selected !"
                self.storage.extra_data['confirm_to_upload'] = _confirm_to_upload          
                return render(request, self.template_name, context)

        # Step 2 - upload
        elif current_step == 'upload': # recheck and save to DB
            form =self.form_list['upload'](request.POST)

            self.upload=True
            
            self.dirname=self.storage.extra_data['dirname'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name  

            # -- Run file processing function
            self.valLog=self.file_process_handler(request, self.dirname, self.filelist, form_data=request.POST, upload=self.upload, appuser=request.user)
            
            # -- Convert valLog results from file processing into HTML table
            if self.valLog.nLogs['Error'] >0 :
                htmlLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns) 
            else:
                htmlLog = self.valLog.get_ashtml(columns=self.html_columns)

            self.storage.extra_data['validation_result'] = htmlLog  
            self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) Uploaded." 

        return self.get_form_step_data(form)

    # Remove uploaded Files
    def _remove_uploaded_files(self):
        _Dirname=self.storage.extra_data['dirname']
        _Filelist = self.storage.extra_data['filelist']
        for f in _Filelist:
            _uploaded_file = os.path.join(_Dirname,f)
            if os.path.isfile(_uploaded_file):
                os.remove(_uploaded_file)

    # Step 3 - upload
    def done(self, form_list, **kwargs):
        
        # Cleanup uploaded Files
        self._remove_uploaded_files()

        # Redirect to the desired page after finishing
        return redirect(self.request.META['HTTP_REFERER'])
    

    def get_context_data(self, form, **kwargs):
        context = super().get_context_data(form=form, **kwargs)

        # save information to context,
        # then display in templates  

        context['step1']=self.name_step1
        current_step = self.steps.current
        context['validation_message'] = self.storage.extra_data.get('validation_message', None)

        if current_step == 'upload_file':
            context['validation_result']="Select Data files to Upload"
        else:
            context['validation_result'] = self.storage.extra_data.get('validation_result', None)
            context['confirm_to_upload']=self.storage.extra_data.get('confirm_to_upload', None)
        #print(f"[ImportHandler_View] {current_step} validation_result: {context['validation_result']}")
        return context