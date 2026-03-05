import os
import shutil
from django import forms
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict

from applib.django.base.views import SuperUserRequiredMixin, WriteUserRequiredMixin
from apputil.utils.files_upload import validate_file, file_location, OverwriteStorage
from applib.logging.validation_log import Validation_Log

# =================================================================
# Utilities Forms
# -----------------------------------------------------------------
class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True

# -----------------------------------------------------------------
class SingleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = False

# class MultipleFileFolderInput(forms.Form):
#     # This field will handle multiple files from the selected directory
#     folder_contents = forms.FileField(
#         widget=forms.ClearableFileInput(attrs={'webkitdirectory': True, 'multiple': True}),
#         required=False # Make it optional if you want
#     )
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

# -----------------------------------------------------------------
# class SingleFolderField(forms.FileField):
    
#     def __init__(self, *args, **kwargs):
#         kwargs.setdefault("widget", SingleFileInput())
#         super().__init__(*args, **kwargs)
       
        
#     def clean(self, data, initial=None):
#         single_file_clean = super().clean
#         if isinstance(data, (list, tuple)):
#             result = [single_file_clean(d, initial) for d in data]
#         else:
#             result = single_file_clean(data, initial)
#         return result

# =================================================================
# Process Step Forms
# -----------------------------------------------------------------
class SelectSingleFile_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    single_file = SingleFileField(label='Select Single File', 
                                  validators=[validate_file], 
                                  required=False,
                                  #help_text="Select a single file"
                                  )
        
    def clean(self):
        cleaned_data = super().clean()
        uploadfiles=[]
        # List of file fields to validate
        file_fields = ['single_file',]
        
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
    multi_files = MultipleFileField(label='Select Multiple Files ', 
                                  validators=[validate_file], 
                                  required=False,
                                  #help_text="Select one or multiple files"
                                  )

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
class SelectSingleFileFolder_StepForm(WriteUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    single_file = SingleFileField(label='Select Single File', 
                                  validators=[validate_file], 
                                  required=False,
                                  #help_text="Select a single file"
                                  )

    multi_files = MultipleFileField(label='Upload Multiple Files', 
                                  validators=[validate_file], 
                                  required=False,
                                  #help_text="Select one or multiple files"
                                  )

    def clean(self):
        cleaned_data = super().clean()
        uploadfiles=[]
        # List of file fields to validate
        file_fields = ['single_file','multi_files']
        
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
        
        print(f" [SelectSingleFileFolder_StepForm] : {file_fields}")
        for e in uploadfiles:
            print(f" [SelectSingleFileFolder_StepForm] UploadFiles: {e}")
        
        return cleaned_data
    
# --------------------------------------------------------------------------------------------------
class Upload_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    upload = forms.BooleanField(initial=False, required=False, help_text="Upload New Data to Database")
    overwrite = forms.BooleanField(initial=False, required=False, help_text="Overwrite Existing Data as well")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #empty_field = models.CharField(max_length=100, blank=True)
        self.fields['upload'].label = "Upload New Data"
        self.fields['overwrite'].label = "Overwrite Existing Data"
        # self.fields['upload'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}
        # self.fields['overwrite'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}

# --------------------------------------------------------------------------------------------------
class Generate_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    generate = forms.BooleanField(initial=False, required=False, help_text="Generate output")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.doc_type = "Report"
        self.fields['generate'].label = "Generate requested Document"

# --------------------------------------------------------------------------------------------------
class Finalize_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    pass

# --------------------------------------------------------------------------------------------------
class Download_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    pass
