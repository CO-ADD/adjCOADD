import os
from django.core.cache import cache

from apputil.utils.validation_log import Validation_Log
from dscreen.models import Screen_Run
from applib.plate.multimode_reader import multimodereader_xls

import logging
logger = logging.getLogger(__name__)



def Upload_ReadOuts_Process(Request, DirName, FileList, RunID=None, upload=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single Vitek PDF, given by:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        RunID: RunID for Screen_Run
        upload : Validation only (False) or Validation and Upload (True)
        appuser : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0

    djRun = Screen_Run.get(RunID)

    valLog = Validation_Log("Upload_ReadOuts")


    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_log('Info', 'Read File', FileList[i]) 
            
            print(f" [Upload_ReadOuts] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            #lstTP = multimodereader_xls(os.path.join(DirName,FileList[i]))
            
            # ,OrgBatchID=OrgBatchID,upload=upload,appuser=appuser,valLog=valLog)

    else:
        print(f"[Upload_ReadOuts] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    
    return(valLog)


def Upload_Motherplates_Process(Request, DirName, FileList, RunID=None, upload=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single Vitek PDF, given by:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        RunID: RunID for Screen_Run
        upload : Validation only (False) or Validation and Upload (True)
        appuser : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0

    djRun = Screen_Run.get(RunID)
    valLog = Validation_Log("Upload_Motherplates")

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_log('Info', 'Read File', FileList[i]) 
            
            print(f" [Upload_Motherplates] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            #lstTP = multimodereader_xls(os.path.join(DirName,FileList[i]))
            
            # ,OrgBatchID=OrgBatchID,upload=upload,appuser=appuser,valLog=valLog)

    else:
        print(f"[Upload_Motherplates] No Xlsx to process in {DirName}  ")

    valLog.select_unique()

    #    if upload:
    #         update_screenrun_summary(djRun)
    #         djRun.save(**kwargs)

    return(valLog)
