import os
from django.core.cache import cache

#from apputil.utils.validation_log import Validation_Log
from applib.logging.validation_log import Validation_Log
from dscreen.models import Screen_Run
from applib.plate.multimode_reader import multimodereader_xls
from applib.plate.masterplates import read_motherplate_prepsheet_xls
from dscreen.utils.summary import update_screenrun_summary

import logging
logger = logging.getLogger(__name__)


#-----------------------------------------------------------------------------------
def Summary_ScreenRun_Process(Request, RunID, upload=False, overwrite=False, appuser=None):
    djRun = Screen_Run.get(RunID)
    print(f" [ScreenRun] Update Summary [{djRun}]")
    update_screenrun_summary(djRun)
    djRun.save()


def Upload_ReadOuts_Process(Request, DirName, FileList, RunID=None, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single File:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        RunID: RunID for Screen_Run
        upload : Validation only (False) or Validation and Upload (True)
        overwrite : On Upload overwrite existing data
        appuser : User Instance of user uploading

    """
    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0
    nUploads = 0

    djRun = Screen_Run.get(RunID)
    valLog = Validation_Log("Upload_ReadOuts")

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_info('Read Readout File', FileList[i]) 
            
            print(f" [Upload_ReadOuts] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            lstTP = multimodereader_xls(os.path.join(DirName,FileList[i]),valLog=valLog)
            
            for _tp in lstTP:
                print(f" [Upload_ReadOuts] Validating: {_tp['plate']} ")
 
                validStatus = True
                validDict = {}
                _tp['plate'].run_id = djRun
                _tp['plate'].set_defaults_model()
                validDict = _tp['plate'].validate_model(WellData=False, verbose = 0)
                if validDict:
                    validStatus = False
                    for c in validDict:
                        print(f" [Upload_ReadOuts] validDict: {c} ")
                        #valLog.add_error('')
                        
                if upload and validStatus:
                    if _tp['new'] or overwrite:
                        _tp['plate'].save(verbose=0)
                        print(f" [Upload_ReadOuts] Saving: {_tp['plate_id']} [Overwrite: {overwrite}]")
                        nUploads += 1
        if upload:
            if len(lstTP)-nUploads > 0:                
                valLog.add_warning("Partial Upload",
                                f"RunID: {djRun.run_id}", 
                                f"Testplates {nUploads} of {len(lstTP)}",
                                "")
            else:
               valLog.add_info("Successful Upload",
                                f"RunID: {djRun.run_id}", 
                                f"Testplates {nUploads} of {len(lstTP)}",
                                "")
                
    else:
        print(f"[Upload_ReadOuts] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    
    return(valLog)


def Upload_Motherplates_Process(Request, DirName, FileList, RunID=None, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single File:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        RunID: RunID for Screen_Run
        upload : Validation only (False) or Validation and Upload (True)
        overwrite : On Upload overwrite existing data
        appuser : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0
    nUploads = 0
    
    djRun = Screen_Run.get(RunID)
    valLog = Validation_Log("Upload_Motherplates")

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_info('Read PlatePrep File', FileList[i],"[MotherPlates]") 
            
            print(f" [Upload_Motherplates] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            lstMP = read_motherplate_prepsheet_xls(os.path.join(DirName,FileList[i]),valLog=valLog)

            for _mp in lstMP:
                print(f" [Upload_Motherplates] Validating: {_mp['plate']} ")
                validStatus = True
                validDict = {}
                
                _mp['plate'].run_id = djRun
                _mp['plate'].set_defaults_model()

                validDict = _mp['plate'].validate_model(WellData=True, verbose = 0)
                if validDict:
                    validStatus = False
                    for c in validDict:
                        print(f" [Upload_MotherPlates] validDict: {c} ")

                if upload and validStatus:
                    if _mp['new'] or overwrite:
                        _mp['plate'].save(verbose=0)
                        print(f" [Upload_MotherPlates] Saving: {_mp['plate']} [Overwrite: {overwrite}]")
                        nUploads += 1
        if upload:
            if len(lstMP)-nUploads > 0:                
                valLog.add_warning("Partial Upload",
                                f"RunID: {djRun.run_id}", 
                                f"Testplates {nUploads} of {len(lstMP)}",
                                "")
            else:
               valLog.add_info("Successful Upload",
                                f"RunID: {djRun.run_id}", 
                                f"Testplates {nUploads} of {len(lstMP)}",
                                "")
    else:
        print(f"[Upload_Motherplates] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    return(valLog)

def Upload_TestplateList_Process(Request, DirName, FileList, RunID=None, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single File:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        FileList : XlsName without FolderName
        RunID: RunID for Screen_Run
        upload : Validation only (False) or Validation and Upload (True)
        overwrite : On Upload overwrite existing data
        appuser : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0
    nUploads = 0
    
    djRun = Screen_Run.get(RunID)
    valLog = Validation_Log("Upload_TestplateList")

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_info('Read PlatePrep File', FileList[i],"[TestPlateList]") 
            
            print(f" [Upload_TestplateList] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            lstMP = read_motherplate_prepsheet_xls(os.path.join(DirName,FileList[i]),valLog=valLog)

        #     for _mp in lstMP:
        #         print(f" [Upload_MothUpload_Testplateserplates] Validating: {_mp['plate']} ")
        #         validStatus = True
        #         validDict = {}
                
        #         _mp['plate'].run_id = djRun
        #         _mp['plate'].set_defaults_model()

        #         validDict = _mp['plate'].validate_model(WellData=True, verbose = 0)
        #         if validDict:
        #             validStatus = False
        #             for c in validDict:
        #                 print(f" [Upload_Testplates] validDict: {c} ")

        #         if upload and validStatus:
        #             if _mp['new'] or overwrite:
        #                 _mp['plate'].save(verbose=0)
        #                 print(f" [Upload_Testplates] Saving: {_mp['plate']} [Overwrite: {overwrite}]")
        #                 nUploads += 1
        # if upload:
        #     if len(lstMP)-nUploads > 0:                
        #         valLog.add_warning("Partial Upload",
        #                         f"RunID: {djRun.run_id}", 
        #                         f"Testplates {nUploads} of {len(lstMP)}",
        #                         "")
        #     else:
        #        valLog.add_info("Successful Upload",
        #                         f"RunID: {djRun.run_id}", 
        #                         f"Testplates {nUploads} of {len(lstMP)}",
                                # "")
    else:
        print(f"[Upload_TestplateList] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    return(valLog)
