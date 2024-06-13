'''
View for uploading Vitek PDFs
'''
import os
import pandas as pd

from django import forms
from apputil.utils.form_wizard_tools import ImportHandler_View, SelectSingleFile_StepForm, Upload_StepForm, Finalize_StepForm
from dorganism.models import Taxonomy
from apputil.models import ApplicationUser, Dictionary
from apputil.utils.validation_log import Validation_Log
from dcell.models import Cell, Cell_Batch
from dcell.utils.upload_cell import upload_Cells_fromXls, Cell_fromDict
# from ddrug.utils.vitek import upload_VitekPDF_Process

import logging
logger = logging.getLogger(__name__)

#=================================================================================================
# Vitek Import Wizard
#=================================================================================================

# -----------------------------------------------------------------
# class VitekValidation_StepForm(Upload_StepForm):
# # -----------------------------------------------------------------
#     orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',
#                                        widget=forms.Select(attrs={'class': 'form-control'}), 
#                                        required=False, 
#                                        help_text='(**Optional, if choose an Organism Batch ID, it will be applied in all uploaded files)',
#                                        queryset=Organism_Batch.objects.filter(astatus__lte=0), )
#     field_order = ['orgbatch_id', 'confirm']

# # -----------------------------------------------------------------
# class Import_VitekView(ImportHandler_View): 
# # -----------------------------------------------------------------    
#     name_step1="Upload"
#     form_list = [
#         ('select_file', SelectFile_StepForm),
#         ('upload', VitekValidation_StepForm),
#         ('finalize', Finalize_StepForm),
#     ]
#     template_name = 'ddrug/importhandler_vitek.html'


#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.orgbatch_id=None
    
#     # customize util functions to validate files:
#     # vitek -- upload_VitekPDF_Process
#     def file_process_handler(self, request, *args, **kwargs):
#         try:
#             form_data=kwargs.get('form_data', None)
#         except Exception as err:
#             print(err)
#             return (err)
#         if 'upload-orgbatch_id' in form_data.keys():
#             self.organism_batch=form_data['upload-orgbatch_id'] #get organism_batch  
#             print(self.organism_batch)   
#         valLog=upload_VitekPDF_Process(request, self.dirname, self.filelist, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user) 
#         return(valLog)
# # 

# -----------------------------------------------------------------
# Process View
# -----------------------------------------------------------------

#=================================================================================================
# Cell Import Wizard
#=================================================================================================

class Import_CellView(ImportHandler_View):
    
    name_step1="Upload"
    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dcell/cell/importhandler_cell.html'
    
    # customize util functions to validate files:
    # vitek -- upload_VitekPDF_Process
    def file_process_handler(self, request, *args, **kwargs):
        try:
            form_data=kwargs.get('form_data', None)
        except Exception as err:
            print(f"[Import_CellView] {err}")
        
        valLog=upload_Cells_Process(request, self.dirname, self.filelist, upload=self.upload, appuser=request.user) 
        return(valLog)
        
# ---------------------------------------------------------------------------------
def upload_Cells_Process(Request, DirName, FileList, upload=False,appuser=None):
    """
    Uploads (upload=True) the data from a single Excel File, given by:
        Request  : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName  : FolderName
        FileList : Excel Workbook name without FolderName
        upload   : Validation only (False) or Validation and Upload (True)
        appuser  : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0

    # Get Xlsx Files in case none given
    if nFiles == 0:
        if DirName:
            DirFiles = os.listdir(DirName)
            FileList = [f for f in DirFiles if f.endswith(".xlsx")]
            nFiles = len(FileList)

    valLog = Validation_Log("upload_Cell_Xls")

    if nFiles > 0:
        for i in range(nFiles):
            logger.info(f"[upload_Cell_Xls_Process] {i+1:3d}/{nFiles:3d} - {FileList[i]}   [{appuser}] ")
            Cell_fromDict(DirName,FileList[i],upload=upload,appuser=appuser,valLog=valLog)

    else:
        logger.info(f"[upload_Cell_Xlsx] NO Excel files to process in {DirName}  ")

    valLog.select_unique()
    
    return(valLog)

    pass

# ---------------------------------------------------------------------------------
# util function to parse and validation drug excel files
# def upload_Cells_fromXls(DirName, FileName, upload=False,appuser=None,valLog=None):

#     if not valLog:
#         valLog = Validation_Log(FileName)

#     print(f"[imp_Cell_fromXlsx]: {DirName} {FileName}")
#     df = pd.read_excel(os.path.join(DirName,FileName))
#     df.columns = [c.lower() for c in df.columns]
#     print(f"[upload_Cell_Xls] {df.columns}")

#     dict = df.to_dict('records')

#     for c in dict:
#         djCell = upload_Cell_fromDict(c,valLog)

    #     e['organism_name'] = Taxonomy.get(e['organism_name'])
    #     e['biologist'] = ApplicationUser.get(e['biologist'])
    #     e['mta_status'] = Dictionary.get(e['mta_status'])
    #     e.pop('old_id')
    #     cell = Cell(**e)
    #     cell.save()
    #     print(cell.cell_names)
    
    #print(dict)
#-----------------------------------------------------------------------------------
   