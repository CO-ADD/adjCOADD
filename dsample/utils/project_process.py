import os
import pandas as pd
from decimal import Decimal
from django.core.cache import cache

from applib.logging.validation_log import Validation_Log
from dsample.models import Project
from dcollab.models import Collab_Group, Collab_User, Collab_Membership, Organisation
from dplate.models import MasterPlate, MasterWell


from applib.plate.stockprep import  read_Stock_Prepsheet_XLS
from applib.project.import_project import get_CompoundSubmisssion_xlsx, parse_SampleInfo_Sheet, parse_ContactInfo_Sheet
from dsample.utils.summary import update_project_summary
from dsample.utils.compounds import Upload_COADD_Compound

from django.conf import settings
import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------
def Load_Project_Process(Request, DirName, FileList,
                         ProjectID=None, UploadContent=['Samples','Contacts'], 
                         upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0
    nUploads = 0

    valLog = Validation_Log("Upload_Project")

    if nFiles > 0:
        for i in range(nFiles):
            
            if settings.DEBUG:
                print(f" [Upload_Project] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{appuser}] ")
                
            _dictSheets = get_CompoundSubmisssion_xlsx(os.path.join(DirName,FileList[i]), valLog=valLog)

            # - Project ----------------------------------------------
            djPrj = None
            if ProjectID:
                djPrj = Project.get(ProjectID)
                if djPrj is None:
                    djPrj = Project()
                    djPrj.project_id = ProjectID
            else:
                djPrj = Project()

            # - Samples ----------------------------------------------
            if 'Samples' in UploadContent and _dictSheets['Samples'] is not None:
                _Samples = parse_SampleInfo_Sheet(_dictSheets['Samples'],valLog=valLog)
                for _smp in _Samples:
                    _smp['project_id'] = str(djPrj)
                    Upload_COADD_Compound(_smp, upload=upload, overwrite=overwrite, valLog=valLog)
                
            # - Contacts / Project Membership ------------------------
            if 'Contacts' in UploadContent and _dictSheets['Contacts'] is not None:
                _Contacts,_PrjTitle = parse_ContactInfo_Sheet(_dictSheets['Contacts'],valLog=valLog)
                
                for key in _Contacts:
                    djUsr = Collab_User.get(None,_Contacts[key]['email'])
                    djOrg = Organisation.get_bysimilarity(_Contacts[key]['organisation'])
                    
                    if djOrg is None:
                        djOrg = Organisation()
                        djOrg.organisation_name = _Contacts[key]['organisation']
                        valLog.add_warning("New Organisation",_Contacts[key]['organisation'])
                    
                    if djUsr is None:
                        djUsr = Collab_User()   
                        djUsr.email = _Contacts[key]['email'] 
                        djUsr.first_name = _Contacts[key]['first_name'] 
                        djUsr.last_name = _Contacts[key]['last_name']
                        djUsr.organisation_id = djOrg 
                        valLog.add_warning("New Collaborator",f"{djUsr}")
                    
                    if 'PI' == key:
    
                        djGrp = Collab_Group.get(None,Code=None, PI_ID=djUsr.user_id)
                        if djGrp is None:
                            djGrp = Collab_Group()
                            djGrp.group_code = f"{_Contacts[key]['first_name'][0]}{_Contacts[key]['last_name']}{djOrg.organisation_code}"
                            djGrp.organisation_id = djOrg
                            valLog.add_warning("New Group",f"{djGrp}")
                              
                        print(f" [Upload_Project] {djGrp}")            
                    print(f" [Upload_Project] {djUsr} [{key}] {djOrg}")            
                        
            valLog.show()
             
                
    else:
        print(f" [Upload_Project] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    return(valLog)

#-----------------------------------------------------------------------------------
def Summary_Project_Process(Request, ProjectID, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    djPrj = Project.get(ProjectID)
    print(f" [Project] Update Summary [{djPrj}]")
    update_project_summary(djPrj)
    djPrj.save()

#-----------------------------------------------------------------------------------
def Upload_StockPrep_Process(Request, DirName, FileList, ProjectID=None, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single File:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        ProjectID: ProjectID for Project
        upload : Validation only (False) or Validation and Upload (True)
        overwrite : On Upload overwrite existing data
        appuser : User Instance of user uploading
    """
    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0
    nUploads = 0

    djPrj = Project.get(ProjectID)
    valLog = Validation_Log("Upload_StockPrep")

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_info('Read StockPrep File', FileList[i]) 

            if settings.DEBUG:
                print(f" [Upload_StockPrep] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djPrj}]  [{appuser}] ")
            lstMP = read_Stock_Prepsheet_XLS(os.path.join(DirName,FileList[i]),valLog=valLog)

            # Process MP's
            for _mpid in lstMP:
                print(f" [Upload_StockPrep] Validating: {lstMP[_mpid]['plate']} ")
                validStatus = True
                validDict = {}
                
                lstMP[_mpid]['plate'].set_defaults_model()

                validDict = lstMP[_mpid]['plate'].validate_model(WellData=True, verbose = 0)
                if validDict:
                    validStatus = False
                    for c in validDict:
                        print(f" [Upload_StockPrep] validDict: {c} ")

                if upload and validStatus:
                    if lstMP[_mpid]['new'] or overwrite:
                        print(f" [Upload_StockPrep] Saving: {lstMP[_mpid]['plate']} [Overwrite: {overwrite}]")
                        # for w in lstMP[_mpid]['plate'].wells:
                        #     _w = lstMP[_mpid]['plate'].wells[w]
                        #     print(f" {_w} {_w.barcode} {_w.cmpbatch_id}")
                        lstMP[_mpid]['plate'].save(verbose=0)
                        nUploads += 1
        if upload:
            if len(lstMP)-nUploads > 0:                
                valLog.add_warning("Partial Upload",
                                f"ProjectID: {djPrj.project_id}", 
                                f"StockPlates {nUploads} of {len(lstMP)}",
                                "")
            else:
               valLog.add_info("Successful Upload",
                                f"ProjectID: {djPrj.project_id}", 
                                f"StockPlates {nUploads} of {len(lstMP)}",
                                "")

    else:
        print(f"[Upload_StockPrep] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    return(valLog)




