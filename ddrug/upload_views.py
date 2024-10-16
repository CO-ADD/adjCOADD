'''
View for uploading Vitek PDFs
'''
import pandas as pd

from django import forms
from apputil.utils.form_wizard_tools import ImportHandler_View, SelectMultipleFiles_StepForm, Upload_StepForm, Finalize_StepForm
from dorganism.models import Organism_Batch
from ddrug.utils.vitek import upload_VitekPDF_Process

#=================================================================================================
# Vitek Import Wizard
#=================================================================================================

# -----------------------------------------------------------------
class VitekValidation_StepForm(Upload_StepForm):
# -----------------------------------------------------------------
    orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',
                                       widget=forms.Select(attrs={'class': 'form-control'}), 
                                       required=False, 
                                       help_text='(**Optional, if choose an Organism Batch ID, it will be applied in all uploaded files)',
                                       queryset=Organism_Batch.objects.filter(astatus__lte=0), )
    field_order = ['orgbatch_id', 'confirm']

# -----------------------------------------------------------------
class Import_VitekView(ImportHandler_View):
# -----------------------------------------------------------------    
    name_step1="Upload"
    form_list = [
        ('select_file', SelectMultipleFiles_StepForm),
        ('upload', VitekValidation_StepForm),
        ('finalize', Finalize_StepForm),
    ]
    template_name = 'ddrug/importhandler_vitek.html'


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.orgbatch_id=None
    
    # customize util functions to validate files:
    # vitek -- upload_VitekPDF_Process
    def file_process_handler(self, request, *args, **kwargs):
        try:
            form_data=kwargs.get('form_data', None)
        except Exception as err:
            print(err)
            return (err)
        if 'upload-orgbatch_id' in form_data.keys():
            self.organism_batch=form_data['upload-orgbatch_id'] #get organism_batch  
            print(self.organism_batch)   
        valLog=upload_VitekPDF_Process(request, self.dirname, self.filelist, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user) 
        return(valLog)
# 

# -----------------------------------------------------------------
# Process View
# -----------------------------------------------------------------


#=================================================================================================
# Drug Import Wizard
#=================================================================================================

class Import_DrugView(ImportHandler_View):
    
    name_step1="Upload"
    form_list = [
        ('select_file', SelectMultipleFiles_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]
    template_name = 'ddrug/importhandler_drug.html'
    
    # customize util functions to validate files:
    # vitek -- upload_VitekPDF_Process
    def file_process_handler(self, request, *args, **kwargs):
        try:
            form_data=kwargs.get('form_data', None)
        except Exception as err:
            print(err)
        
        valLog=imp_Drug_fromXlsx(request, self.dirname, self.filelist, upload=self.upload, appuser=request.user) 
        return(valLog)
        
        
# ---------------------------------------------------------------------------------
# util function to parse and validation drug excel files
def imp_Drug_fromXlsx(request, *args, **kwargs):
    pass
#-----------------------------------------------------------------------------------
   