import os
import pandas as pd
from decimal import Decimal
from django.core.cache import cache

from applib.logging.validation_log import Validation_Log
from dsample.models import Project
from dplate.models import MasterPlate, MasterWell


from applib.plate.stockprep import  read_Stock_Prepsheet_XLS
from dsample.utils.summary import update_project_summary

from django.conf import settings
import logging
logger = logging.getLogger(__name__)


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


