import os
from pathlib import Path
import pandas as pd

from apputil.utils.validation_log import Validation_Log
from dorganism.models import Taxonomy
from dcell.models import Cell, Cell_Batch

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import split_StrList, append_StrList

import logging
logger = logging.getLogger(__name__)

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
# ----------------------------------------------------------------------------------------------------

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

    valLog = Validation_Log("upload_Cells_Xls")

    if nFiles > 0:
        for i in range(nFiles):
            logger.info(f"[upload_Cell_Xls_Process] {i+1:3d}/{nFiles:3d} - {FileList[i]}   [{appuser}] ")
            upload_Cells_fromXls(DirName,FileList[i],upload=upload,appuser=appuser,valLog=valLog)

    else:
        logger.info(f"[upload_Cell_Xls_Process] NO Excel files to process in {DirName}  ")

    valLog.select_unique()
    
    return(valLog)

# ----------------------------------------------------------------------------------------------------
def upload_Cells_fromXls(DirName, FileName, upload=False,appuser=None,valLog=None):
    """
    Read Excel Sheet for Cell Uploads 
    """
# ----------------------------------------------------------------------------------------------------
    if not valLog:
        valLog = Validation_Log(FileName)

    #print(f"[upload_Cells_fromXls]: {DirName} {FileName}")

    df = pd.read_excel(os.path.join(DirName,FileName))
    df.columns = [c.lower() for c in df.columns]
    dict = df.to_dict('records')

    for c in dict:
        djCell = Cell_fromDict(c,valLog,for_upload=upload)
        if upload:
            if djCell.VALID_STATUS:
                logger.debug(f" {djCell} <- {FileName}")
                djCell.save(user=appuser)
                c['cell_id'] = str(djCell)
            else:
                valLog.add_log('Error','Cell not validated',f"{c['cell_line']}",'-')    
        print(f" [upload_Cells_fromXls] cell_id: {c['cell_id']}")
        djCellBatch = CellBatch_fromDict(c,valLog,for_upload=upload)
        if upload:
            if djCellBatch.VALID_STATUS:
                logger.debug(f" {djCellBatch} <- {FileName}")
                djCellBatch.save(user=appuser)
            else:
                valLog.add_log('Error','CellBatch not validated',f"{c['cell_line']}",'-')    


# ----------------------------------------------------------------------------------------------------
def Cell_fromDict(iDict,valLog,for_upload=True):
    """
    Create Cell instance from {iDict}
    """
# ----------------------------------------------------------------------------------------------------

    # cell_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Cell ID") 
    # cell_line= models.CharField(max_length=200, blank=True, verbose_name = "Cell Line") 
    # cell_names= models.CharField(max_length=200, blank=True, verbose_name = "Cell Names") 
    # cell_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Cell Notes")
    # cell_code= models.CharField(max_length=30, blank=True, verbose_name = "Cell Code")
    # cell_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    # cell_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    # cell_identification = models.CharField(max_length=512, blank=True, verbose_name = "Cell Identification")
    # cell_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Cell")
    # source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    # source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    # reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    # organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
    #     db_column="organism_name", related_name="%(class)s_organism_name")

    # mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
    #     db_column="mta_status", related_name="%(class)s_mta")
    # mta_notes = models.CharField(max_length=150, blank=True, verbose_name = "MTA Notes")
    # mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    # mta_notes = models.CharField(max_length=512, blank=True, verbose_name = "MTA Notes")

    # collect_tissue = models.CharField(max_length=120, blank=True, verbose_name = "From Tissue/Organ")
    # patient_diagnosis = models.CharField(max_length=120, blank=True, verbose_name = "Patient Diagnosis")
    # patient = models.CharField(max_length=20, blank=True, verbose_name = "Patient Info")

    # biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
    #     db_column="biologist", related_name="%(class)s_biologist")

    # assoc_documents = models.ManyToManyField(Document,verbose_name = "Documents", blank=True,
    #     db_table = "cell_doc", related_name="%(class)s_document")

    # Choice_Dictionary = {
    #     'mta_status':'License_Status',
    #     'cell_type':'Cell_Type',
    #     'cell_panel':'Cell_Panel',
    # }



    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 
    validStatus = True
    vlog_Process = "Upload Cell"

    # -- Remove nan 
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    # -- Find Cell Instance if exist
    djCell = None
    if 'cell_id' in iDict:
        if iDict['cell_id'] is not None:
            djCell = Cell.get(iDict['cell_id'],None)
    if djCell is None:
        djCell = Cell.get(None,iDict['cell_line'])

    if djCell is None:
        djCell = Cell()
        djCell.cell_line = iDict['cell_line']
        valLog.add_log("Info",vlog_Process,f"{iDict['cell_line']}",'New Cell Line','-') 

    # -- Find ForeignKeys
    _OrgName = Taxonomy.get(iDict['organism_name'])
    if _OrgName is None:
        validStatus = False
        valLog.add_log("Error",vlog_Process,f"{iDict['organism_name']} ",'[Organism Name] does not Exists','-')
    else:
        djCell.organism_name = _OrgName

    if 'biologist' in iDict:
        djCell.biologist = ApplicationUser.get(iDict['biologist'])
        if djCell.biologist is None:
            valLog.add_log("Error",vlog_Process,f"{iDict['biologist']} ",'[Biologist]  not found','-')
            validStatus = False

    # -- Special Consideration
    # if 'old_id' in iDict:
    #     if 'cell_notes' in iDict:
    #         iDict['cell_notes'] = append_StrList(iDict['cell_notes'], f"Previous ID: {iDict['old_id']}")
    #     else:
    #         iDict['cell_notes'] = f"[Old_ID]: {iDict['old_id']};"    

    # -- Find Single Dictionaries
    _SngDictFields = ['mta_status']
    for _field in _SngDictFields:
       if _field in iDict:
            if iDict[_field] is not None:
                _dictValue = Dictionary.get(Cell.Choice_Dictionary[_field],iDict[_field],verbose=0)
                if _dictValue is None:
                    valLog.add_log("Error",vlog_Process,f"Value: {iDict[_field]} ",f"[{_field}] not found",'-')
                    validStatus = False
                setattr(djCell, _field, _dictValue)


    # -- Find List of Dictionaries
    _ArrDictFields = ['cell_type','cell_panel']
    for _field in _ArrDictFields:
       if _field in iDict:
            if iDict[_field] is not None:
                _dictList, _errList = Dictionary.get_DictValues_fromStrList(Cell.Choice_Dictionary[_field],iDict[_field],verbose=0)
                if len(_errList) > 0 :
                    for _errValue in _errList:
                        valLog.add_log("Error",vlog_Process,f"Value: {_errValue} ",f"[{_field}] not found",'-')
                        validStatus = False
                setattr(djCell, _field, _dictList)

        
    # -- Assign normal fields
    _AssignFields = ['cell_names','cell_notes','cell_code','cell_identification','cell_origin',
                  'source','source_code','reference',
                  'mta_document',
                  'collect_tissue','patient_diagnosis','patient']
    for _field in _AssignFields:
        if _field in iDict:
            setattr(djCell, _field, iDict[_field])

    # -- Clean and Validate Entry
    print(f"---------------------------------")
    djCell.clean_Fields()
    validStatus = True
    validDict = djCell.validate()
    print(f"{djCell}")
    print(f"{validDict}")

    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log("Warning",vlog_Process,k,validDict[k],'-')

    djCell.VALID_STATUS = validStatus

    return(djCell)

# ----------------------------------------------------------------------------------------------------
def CellBatch_fromDict(iDict,valLog,for_upload=True):
    """
    Create Cell instance from {iDict}
    """
# ----------------------------------------------------------------------------------------------------

    # cellbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "CellBatch ID")
    # cell_id = models.ForeignKey(Cell, null=False, blank=False, verbose_name = "Cell ID", on_delete=models.DO_NOTHING,
    #     db_column="cell_id", related_name="%(class)s_cell_id")
    # previous_batch_id= models.CharField(max_length=20, verbose_name = "Previous CellBatch ID")
    # passage_number= models.CharField(max_length=20, verbose_name = "Passage Number")
    # batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[alphanumeric], verbose_name = "Batch ID")
    # batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")
    # batch_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Quality", on_delete=models.DO_NOTHING,
    #     db_column="batch_quality", related_name="%(class)s_batchquality")
    # quality_source = models.CharField(max_length=150, blank=True, verbose_name = "QC Source")
    # qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC status", on_delete=models.DO_NOTHING,
    #     db_column="qc_status", related_name="%(class)s_qc")
    # qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")
    # stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date") 
    # stock_level = models.CharField(max_length=20, blank=True, verbose_name = "Stock Levels") 
    # biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
    #     db_column="biologist", related_name="%(class)s_biologist")

    iDict =  {k.lower(): v for k, v in iDict.items()} 
    validStatus = True
    vlog_Process = "Upload CellBatch"

    # -- Remove nan 
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    # -- Find CellBatch Instance if exist
    djCellBatch = None
    if 'cellbatch_id' in iDict:
        if iDict['cellbatch_id'] is not None:
            djCellBatch = Cell_Batch.get(iDict['cellbatch_id'])

    if djCellBatch is None:
        djCellBatch = Cell_Batch()
        #djCell.cell_line = iDict['cell_line']
        if for_upload:
            valLog.add_log("Info",vlog_Process,f"{iDict['cell_line']}",'New CellBatch created.','-')
        else: 
            valLog.add_log("Info",vlog_Process,f"{iDict['cell_line']}",'New CellBatch to be created.','-') 

    # -- Find ForeignKeys
    _CellID = None
    if 'cell_id' in iDict:
        if iDict['cell_id'] is not None:
            _CellID = Cell.get(iDict['cell_id'],None)

    if _CellID is None:
        validStatus = False
        if for_upload:
            valLog.add_log("Error",vlog_Process,f"{iDict['cell_id']} ",'[Cell ID] does not exist','-')
        else:
            valLog.add_log("Warning",vlog_Process,f"{iDict['cell_id']} ",'[Cell ID] will be assigned on Upload','-')
    else:
        djCellBatch.cell_id = _CellID

    if 'biologist' in iDict:
        djCellBatch.biologist = ApplicationUser.get(iDict['biologist'])
        if djCellBatch.biologist is None:
            valLog.add_log("Error",vlog_Process,f"{iDict['biologist']} ",'[Biologist]  not found','-')
            validStatus = False

    # -- Find Single Dictionaries
    _SngDictFields = ['qc_status','batch_quality']
    for _field in _SngDictFields:
       if _field in iDict:
            if iDict[_field] is not None:
                _dictValue = Dictionary.get(Cell.Choice_Dictionary[_field],iDict[_field],verbose=0)
                if _dictValue is None:
                    valLog.add_log("Error",vlog_Process,f"Value: {iDict[_field]} ",f"[{_field}] not found",'-')
                    validStatus = False
                setattr(djCellBatch, _field, _dictValue)

        
    # -- Assign normal fields
    _AssignFields = ['passage_number','batch_notes',
                  'quality_source','qc_record','stock_date','stock_level']
    for _field in _AssignFields:
        if _field in iDict:
            setattr(djCellBatch, _field, iDict[_field])

    # -- Clean and Validate Entry
    djCellBatch.clean_Fields()
    validStatus = True
    validDict = djCellBatch.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',vlog_Process,k,validDict[k],'-')

    djCellBatch.VALID_STATUS = validStatus

    return(djCellBatch)