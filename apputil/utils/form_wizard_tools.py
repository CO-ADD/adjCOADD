import os
from django import forms
from django.shortcuts import render
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from django.core.exceptions import ValidationError
from django.utils.datastructures import MultiValueDict
from apputil.utils.views_base import SuperUserRequiredMixin
from apputil.utils.files_upload import validate_file, file_location

class UploadFileForm(SuperUserRequiredMixin, forms.Form):
    multi_files = forms.FileField(label='Select files', widget=forms.ClearableFileInput(attrs={'multiple': True}),validators=[validate_file],  required=False)
    folder_input = forms.FileField(label='Select a folder', widget=forms.ClearableFileInput(attrs={'multiple': True, 'webkitdirectory':True}), validators=[validate_file], required=False)
    
    def clean(self):
        cleaned_data = super().clean()
        print('clean data')

        # List of file fields to validate
        file_fields = ['multi_files','folder_input']
        if not isinstance(self.files, MultiValueDict):
            return cleaned_data

        for field in file_fields:
            files = self.files.getlist(f'upload_file-{field}')
            
            for file in files:
                for validator in self.fields[field].validators:
                    try:
                        validator(file)
                    except ValidationError as e:
                        self.add_error(field, f"{file.name}: {str(e)}")
        return cleaned_data
    

    
class StepForm_1(forms.Form):
    pass

class StepForm_2(forms.Form):
    confirm_to_save = forms.BooleanField(required=True)

class FinalizeForm(forms.Form):
    pass

class ImportHandler_WizardView(SessionWizardView):
    # here add steps name
    step1='step1'
    step2='step2'
    # ...

    form_list = [
        ('upload_file', UploadFileForm),
        # Here adding steps:
        # ('Step1', StepForm_1),
        # ...
        ('finalize', FinalizeForm),
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

        if current_step == 'upload_file':
            # here uploads files and define file name and path
            location = file_location(request)  # define file store path during file process
            files = []
            if 'single_file' in request.FILES:
                files.append(request.FILES['single_file'])
            if 'folder_files' in request.FILES:
                files.extend(request.FILES.getlist('folder_files'))
            # here can add extra functions for files


        elif current_step == 'step1':
            print("2")
            # In this step, you can perform further steps

        elif current_step == 'step2':
            print("3")
            # more steps

        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        print("4")
        # Save data to the database or perform any other final actions
        # Redirect to the desired page after finishing
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        context = super().get_context_data(form=form, **kwargs)
        current_step = self.steps.current

        if current_step == 'step1':
            pass
        #    here can define extra context for step1 result 

        if current_step == 'confirm_validation':
            pass
         #    here can define extra context for step2 result         

        return context
    
        # Use to delete uploaded files
    def delete_file(self, file_name):
        location=file_location(self.request)
        file_full_path=os.path.join(location, file_name)
        print(file_full_path)
        try:
            os.unlink(file_full_path)
            print("removed!")
        except Exception as err:
            raise Exception

    