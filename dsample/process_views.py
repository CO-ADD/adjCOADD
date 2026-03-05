from django import forms

from dsample.models import Project
from dsample.utils.project_process import (Load_Project_Process,
                                           Upload_StockPrep_Process, Summary_Project_Process)
from applib.process.process_forms import Process_View
from applib.process.process_stepforms import SelectSingleFile_StepForm, Finalize_StepForm, Upload_StepForm, Generate_StepForm, SelectSingleFileFolder_StepForm


# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
class UploadSelect_StepForm(Upload_StepForm):
# --------------------------------------------------------------------------------------------------
    contacts = forms.BooleanField(initial=False, required=False, help_text="Upload Contact Information for Project")
    samples = forms.BooleanField(initial=False, required=False, help_text="Upload Compound Information for Project")
    #overwrite = forms.BooleanField(initial=False, required=False, help_text="Overwrite Existing Data as well")
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['contacts'].label = "Contact/Collaborator Information"
        self.fields['samples'].label = "Compound/Sample Information"
        #self.fields['overwrite'].label = "Overwrite Existing Data"
        # self.fields['upload'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}
        # self.fields['overwrite'].error_messages = {'required': 'File(s) contain Errors. Please correct the content of the files'}

# --------------------------------------------------------------------------------------------------
class Load_Project_ProcessView(Process_View):
    process_name = 'Import_Project'
    model = Project
    progress = [0,0]

    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', UploadSelect_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dsample/project_process/load_cmpdsubmission.html'

    select_html  = 'Please select a Compound Submission Excel [xlsx] file<br>'
    select_html += 'Make sure file contains correct <b>Collaborator</b> '
    select_html += 'and <b>Compound</b> informations.'

    upload_html  = 'Please check<br>'
    upload_html += 'Samples - upload Compound Information<br>'
    upload_html += 'Contacts - upload Collaborator Information<br>'

    initial_dict = {
        'select_file': {'instructions':select_html},
        'upload': {'instructions':upload_html},
        'finalize': {'instructions':''},
        }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.project_id=None
        print(f" [{self.process_name}] Init")

    def file_process_handler(self, request, *args, **kwargs):    
        #print(" [Add_TestplateList_ProcessView.file_process_handler]")

        self.samples = False
        self.contacts = False
        self.upload = False
        self.overwrite = False

        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'overwrite' in form_data:
            self.overwrite = form_data['overwrite']

        valLog=Load_Project_Process(request, self.file_dir, self.file_list['single_file'], ProjectID=None,
                                    UploadContent={'Samples':self.samples,'Contacts':self.contacts},
                                    upload=self.upload, overwrite=False,
                                    appuser=request.user)
         
        return(valLog)
    
# --------------------------------------------------------------------------------------------------
class Add_ProjectInfo_ProcessView(Process_View):
    process_name = 'Add_ProjectInfo'
    model = Project
    progress = [0,0]

    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', UploadSelect_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dsample/project_process/load_cmpdsubmission.html'

    select_html  = 'Please select a Compound Submission Excel [xlsx] file<br>'
    select_html += 'Make sure file contains correct <b>Collaborator</b> '
    select_html += 'and <b>Compound</b> informations.'

    upload_html  = 'Please check<br>'
    upload_html += 'Samples - upload Compound Information<br>'
    upload_html += 'Contacts - upload Collaborator Information<br>'

    initial_dict = {
        'select_file': {'instructions':select_html},
        'upload': {'instructions':upload_html},
        'finalize': {'instructions':''},
        }
    

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.project_id=None
        print(f" [{self.process_name}] Init")

    def file_process_handler(self, request, *args, **kwargs):    
        #print(" [Add_TestplateList_ProcessView.file_process_handler]")

        self.samples = True
        self.contacts = True
        self.upload = False
        self.overwrite = False

        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        
        if 'samples' in form_data:
            self.samples = form_data['samples']
        if 'contacts' in form_data:
            self.contacts = form_data['contacts']

        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'overwrite' in form_data:
            self.overwrite = form_data['overwrite']

        valLog=Load_Project_Process(request, self.file_dir, self.file_list['single_file'], ProjectID=self.pk,
                                    UploadContent={'Samples':self.samples,'Contacts':self.contacts},
                                    upload=self.upload, overwrite=self.overwrite,
                                    appuser=request.user)
         
        return(valLog)
    
# --------------------------------------------------------------------------------------------------
class Load_StockPrep_ProcessView(Process_View):
    process_name = 'Upload_StockPrep'
    model = Project
    progress = [0,0]

    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dsample/project_process/load_stockprep.html'

    # select_html  = 'Please select a Excel [xlsx] file from Tecan/BioTek readers'
    # select_html += '\n Make sure file contains correct  <b>TestPlate IDs</b>'

    # upload_html  = 'Please check the TestPlate IDs [<i>Item</i>] for any "New Testplate" [<i>Action</i>]'
    # upload_html += '\n Make sure the IDs are unique and reflect the IDs in <b>TestPLateList</b>'
    # upload_html += '\n In case, correct the IDs in the <b>Readout</b> file and repeat the upload'

    # message_html =[
    #    ('select_file',select_html),
    #    ('upload',upload_html),
    #    ('finalize','') 
    # ]
    def file_process_handler(self, request, *args, **kwargs):    
        #print(" [Add_TestplateList_ProcessView.file_process_handler]")

        self.upload = False
        self.overwrite = False
                
        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'overwrite' in form_data:
            self.overwrite = form_data['overwrite']

        valLog=Upload_StockPrep_Process(request, self.file_dir, self.file_list['single_file'], ProjectID=self.pk, 
                                           upload=self.upload, overwrite=self.overwrite,
                                           appuser=request.user) 
 
        return(valLog)

    def file_process_finalizer(self, request, pk):
        Summary_Project_Process(request, pk)
