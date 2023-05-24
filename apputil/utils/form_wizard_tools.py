'''
Steps process driven by form submission requests. for imporing data with PDFs
-uploading PDFs
-parsing, validating and saving to DB 
'''
import os
from django import forms
from django.shortcuts import render
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict
from apputil.utils.views_base import SuperUserRequiredMixin
from apputil.utils.files_upload import validate_file,file_location


# --------------------------------------------------------------------------------------------------
class SelectFile_StepForm(SuperUserRequiredMixin, forms.Form):
# --------------------------------------------------------------------------------------------------
    multi_files = forms.FileField(label='Select one or multiple files', 
                                  widget=forms.ClearableFileInput(attrs={'multiple': True,}),
                                  validators=[validate_file], 
                                  required=False)

    def clean(self):
        cleaned_data = super().clean()
        print('clean data')
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
class StepForm_1(forms.Form):
# --------------------------------------------------------------------------------------------------
    confirm = forms.BooleanField(required=True, help_text="Confirm to upload Data")
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['confirm'].error_messages = {'required': 'Contains Error, please correct firstly'}


# class StepForm_2(forms.Form):

# --------------------------------------------------------------------------------------------------
class Finalize_StepForm(forms.Form):
# --------------------------------------------------------------------------------------------------
    pass
    #log_entry = forms.CharField(widget=forms.Textarea, required=False)

# --------------------------------------------------------------------------------------------------
class ImportHandler_WizardView(SuperUserRequiredMixin,SessionWizardView):
# --------------------------------------------------------------------------------------------------
    # here add steps name
    step1='upload'
    # step2='step2'
    # ...

    form_list = [
        ('select_file', SelectFile_StepForm),
        ('upload', Upload_StepForm),
        # Here adding steps:
        # ('Step1', StepForm_1),
        # ...
        ('finalize', Finalize_StepForm),
    ]

    template_name = None

    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Here add process variables
    # 
    # process_step should be customized per application data
    def process_step(self, form):
        # Here control steps
        current_step = self.steps.current
        request = self.request

        if current_step == 'select_file':
            # here uploads files and define file name and path
            location = file_location(request)  # define file store path during file process
            files = []
            if 'single_file' in request.FILES:
                files.append(request.FILES['single_file'])
            if 'folder_files' in request.FILES:
                files.extend(request.FILES.getlist('folder_files'))
            # here can add extra functions for files


        elif current_step == 'upload':
            print("Upload")
            # In this step, you can perform further steps

        #     # more steps

        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
      
        # Save data to the database or perform any other final actions
        # Redirect to the desired page after finishing
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        context = super().get_context_data(form=form, **kwargs)
        current_step = self.steps.current

        if current_step == 'upload':
            pass
        #    here can define extra context for upload result 
     

        return context
    
        # Use to delete uploaded files
    def delete_file(self, file_name):
        location=file_location(instance=self.request.user)
        file_full_path=os.path.join(location, file_name)
        print(file_full_path)
        try:
            os.unlink(file_full_path)
            print("removed!")
        except FileNotFoundError:
            print(f"File {file_path} does not exist.")
        except Exception as err:
            raise Exception

    
