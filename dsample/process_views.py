
from dsample.models import Project
from dsample.utils.project_process import (Load_Project_Process,
                                           Upload_StockPrep_Process, Summary_Project_Process)
from applib.process.process_forms import Process_View
from applib.process.process_stepforms import SelectSingleFile_StepForm, Finalize_StepForm, Upload_StepForm, Generate_StepForm, SelectSingleFileFolder_StepForm


# --------------------------------------------------------------------------------------------------
class Load_Project_ProcessView(Process_View):
    process_name = 'Import_Project'
    model = Project
    progress = [0,0]

    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dsample/project_process/load_cmpdsubmission.html'

    select_html  = 'Please select a Compound Submission Excel [xlsx] file '
    select_html += '\n Make sure file contains correct <b>Collaborator</b> '
    select_html += '\n and <b>Compound</b> informations.'

    upload_html  = 'Please check the TestPlate IDs [<i>Item</i>] for any "New Testplate" [<i>Action</i>]'
    upload_html += '\n Make sure the IDs are unique and reflect the IDs in <b>TestPLateList</b>'
    upload_html += '\n In case, correct the IDs in the <b>Readout</b> file and repeat the upload'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.project_id=None
        print(f" [{self.process_name}] Init")

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

        valLog=Load_Project_Process(request, self.file_dir, self.file_list['single_file'],
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
