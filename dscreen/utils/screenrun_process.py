import os
from django.core.cache import cache

#from apputil.utils.validation_log import Validation_Log
from applib.logging.validation_log import Validation_Log
from dscreen.models import Screen_Run
from dplate.models import Labware, TestPlate, TestWell, MasterPlate

from applib.plate.multimode_reader import multimodereader_xls
#from applib.plate.masterplates import read_motherplate_prepsheet_xls
from applib.plate.plateprep import read_Motherplates_Prepsheet_XLS, read_TestPlateList_Prepsheet_XLS, get_BarcodeScans, gen_Motherplates_PSPrep
from applib.plate.testplates import add_mother_to_testplate
from applib.bio.doseresponse import process_testplate_doseresponse
from dscreen.utils.summary import update_screenrun_summary
from adjcoadd.constants import DR_CLASSES

from django.conf import settings
import logging
logger = logging.getLogger(__name__)


#-----------------------------------------------------------------------------------
def Summary_ScreenRun_Process(Request, RunID, upload=False, overwrite=False, appuser=None):
#-----------------------------------------------------------------------------------
    djRun = Screen_Run.get(RunID)
    print(f" [ScreenRun] Update Summary [{djRun}]")
    update_screenrun_summary(djRun)
    djRun.save()


#-----------------------------------------------------------------------------------
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
            
            if settings.DEBUG:
                print(f" [Upload_ReadOuts] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
            lstTP = multimodereader_xls(os.path.join(DirName,FileList[i]),valLog=valLog)
            
            for _tp in lstTP:
                if settings.DEBUG:
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

#-----------------------------------------------------------------------------------
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
            lstMP = read_Motherplates_Prepsheet_XLS(os.path.join(DirName,FileList[i]),valLog=valLog)

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

#-----------------------------------------------------------------------------------
def Upload_TestplateList_Process(Request, DirName, FileList, RunID=None, 
                                 upload=False, apply_mp=False, overwrite=False, appuser=None):
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


    logNumbers = {'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                  'Valid Plates':0, 'Rejected Plates':0, 'Failed Plates':0,
                  'Processed AssayData':0,'Uploaded AssayData':0,
                  'Inhibition AssayData':0, 'MIC AssayData':0, 'CC50 AssayData':0, 'HC50 AssayData':0, 
                  'Empty':0}

    if nFiles > 0:
        for i in range(nFiles):
            valLog.add_info('Read PlatePrep File', FileList[i],"[TestPlateList, Assay]") 
            if settings.DEBUG:
                print(f" [Upload_TestplateList] {i+1:3d}/{nFiles:3d} - {FileList[i]}  [{djRun}]  [{appuser}] ")
                
            lstTP,lstAss = read_TestPlateList_Prepsheet_XLS(os.path.join(DirName,FileList[i]),apply_mp=apply_mp,valLog=valLog)

            if valLog.if_noerrors():
                _MPS = {}
                for _tpid in lstTP:
                    _tp = lstTP[_tpid]['plate']
                    _new = lstTP[_tpid]['plate']
                    validStatus = True
                    if settings.DEBUG:
                        print(f" [Upload_TestplateList] {_tp.plate_id} Validating:  {_tp.result_type} {_tp.assay_id} ")
                    validDict = _tp.validate_model(WellData=False, verbose = 0)
                    if validDict:
                        validStatus = False
                        for c in validDict:
                            print(f" [Upload_TestPlateList] validDict: {c} ")

                    if apply_mp:
                        # Apply MP -> Add Compounds
                        # -------------------------
                        if hasattr(_tp,'motherplate_ids'):
                            if settings.DEBUG:
                                print(f" [Upload_TestplateList] {_tp.plate_id} Apply MP: {_tp.motherplate_ids}")

                            _tp.clear_cmpbatch_data()
                            mp_ids = getattr(_tp,'motherplate_ids')
                            for mp in mp_ids:
                                if mp != '-':
                                    if mp not in _MPS:
                                        _MPS[mp] = MasterPlate.get(mp,WellData=True)
                                    add_mother_to_testplate(_MPS[mp],_tp)
                            
                            _tp.update_n('n_samples')
                        else:
                            validStatus = False
                            valLog.add_error("No MIssing MotherPlates",_tpid,"No MotherPlate_IDs","Correct TestPlateList")

                        # Analyze Testplate 
                        # -----------------
                        if validStatus and _tp.n_wells > 0 and _tp.n_reads > 0 and _tp.control_layout:
                            #print(f"{_tp.plate_id} {_tp.n_wells} {_tp.control_layout}")
        
                            if _tp.apply_layout() > 0:
                                # Analyze Testplate -> Inhibition
                                # -------------------------------
                                _tp.calc_inhibition()
                                if settings.DEBUG:
                                    print(f" [Upload_TestPlateList] {_tpid} {_tp.plate_quality} {_tp.zfactor} {_tp.n_inhibitions}")

                                _dr_list = []
                                if str(_tp.plate_quality) == 'Valid':
                                    logNumbers['Valid Plates'] += 1
                                    if str(_tp.result_type) in DR_CLASSES:
                                        # Analyze Testplate -> DR AssayData                                            
                                        # ---------------------------------
                                        _dr_list = process_testplate_doseresponse(_tp)
                                                                                                    
                                        # Analyze Testplate -> AssayData
                                        # ------------------------------             
                                        logNumbers['Processed AssayData'] += len(_dr_list)
                                        for _dr in _dr_list:
                                            logNumbers[f'{_dr.dr_type} AssayData'] += 1

                                        if settings.DEBUG:
                                            print(f" [Upload_TestplateList] {_tpid} DR: {_tp.result_type} {len(_dr_list)}")
                                            
                                    elif str(_tp.result_type) == 'Inhibition':
                                        logNumbers[f'Inhibition AssayData'] += _tp.n_inhibitions

                                        if settings.DEBUG:
                                            print(f" [Upload_TestplateList] {_tpid} SC: {_tp.result_type} {_tp.n_inhibitions}")
                                    else:
                                        if settings.DEBUG:
                                            print(f" [Upload_TestplateList] {_tpid} ??: {_tp.result_type}")
                                        

                                elif str(_tp.plate_quality) == 'Rejected':
                                    logNumbers['Rejected Plates'] += 1
                                else:
                                    logNumbers['Failed Plates'] += 1
                                
                                # print(f"PosCntr: {_tp.poscontrol_stats}")
                                # print(f"NegCntr: {_tp.negcontrol_stats}")
                                # print(f"Sample: {_tp.sample_stats}")
                                # print(f"Edge: {_tp.edge_stats}")                                
                                # validDict= _tp.validate_fields()
                                # print(validDict)

 
                        else:
                            validStatus = False
                            valLog.add_error("Calc Error",_tp.plate_id,"Either no Wells, Readout or Layout","Correct TestPlateList")

                    if upload and validStatus:
                        if settings.DEBUG:
                            print(f" [Upload_TestPlateList] {_tpid} Updating: {_tp.plate_quality} {_tp.zfactor} {_tp.n_inhibitions}")
                        _tp.save()
                        for _dr in _dr_list:
                            _dr.save_assaydata(overwrite=True)
                        logNumbers['Uploaded Plates'] += 1
                        logNumbers['Uploaded AssayData'] += len(_dr_list)

                # Message valLog
                if apply_mp:

                    #TestPlates
                    valLog.add_info('TestPlates',f"Valid: {logNumbers['Valid Plates']} ","Testplates pass QC","Confirm Testplate")
                    if logNumbers['Rejected Plates'] > 0:
                        valLog.add_warning('TestPlates',f"Rejected: {logNumbers['Rejected Plates']} ","Testplates Rejected by ZFactor QC","Check Testplate/Layout")
                    else:
                        valLog.add_info('TestPlates',f"Rejected: {logNumbers['Rejected Plates']} ","Testplates Rejected by ZFactor QC")

                    if logNumbers['Failed Plates'] > 0:
                        valLog.add_warning('TestPlates',f"Failed: {logNumbers['Failed Plates']} ","Testplates fail ","Check Testplate")
                    else:
                        valLog.add_info('TestPlates',f"Failed: {logNumbers['Failed Plates']} ","Testplates fail ")

                    # AssayData
                    valLog.add_info('SC AssayData',f"Inhibition: {logNumbers['Inhibition AssayData']} ","Inhibition AssayData generated ")
                    valLog.add_info('DR AssayData',f"Processed: {logNumbers['Processed AssayData']} ","Dosereponse AssayData generated ")
                    valLog.add_info('DR AssayData',f"MIC: {logNumbers['MIC AssayData']} ","MIC AssayData generated ")
                    valLog.add_info('DR AssayData',f"CC50: {logNumbers['CC50 AssayData']} ","CC50 AssayData generated ")
                    valLog.add_info('DR AssayData',f"HC50: {logNumbers['HC50 AssayData']} ","HC50 AssayData generated ")
                    if logNumbers['Processed AssayData'] > 0:
                        valLog.add_warning('TestPlates',f"Failed: {logNumbers['Failed Plates']} ","Testplates fail ","Check Testplate")
                    else:
                        valLog.add_info('TestPlates',f"Failed: {logNumbers['Failed Plates']} ","Testplates fail ")

                else:
                    valLog.add_warning('TestPlates',f"No Analysis",f"No Compound, Layout or Analysis applied","Select Apply_MP")

                if upload:
                    valLog.add_info('TestPlates',f"Uploaded: {logNumbers['Uploaded Plates']} ","Testplates uploaded")
                    _not_uploaded = logNumbers['Processed Plates']- logNumbers['Uploaded Plates']
                    if  _not_uploaded > 0:
                        valLog.add_warning('TestPlates',f"Not Uploaded: {_not_uploaded} ","Testplates not uploaded","Check Testplate")
                    valLog.add_info('DR AssayData',f"Uploaded: {logNumbers['Uploaded AssayData']} ","Dosereponse AssayData uploaded")

                else:
                    valLog.add_warning('TestPlates',f"Uploaded: 0",f"No Testplate uploaded to database","Select Upload")

    else:
        print(f"[Upload_TestplateList] No Xlsx to process in {DirName}  ")

    valLog.select_unique()
    return(valLog)

#-----------------------------------------------------------------------------------
def Gen_Masterplates_Process(Request, DirName, PrepFileList, RackFileList, RunID=None, 
                                 generate=False,  appuser=None):
#-----------------------------------------------------------------------------------
    
    if PrepFileList:
        nFiles = len(PrepFileList)
    else:
        nFiles = 0

    if RackFileList:
        nRacks = len(RackFileList)
    else:
        nRacks = 0
    
    djRun = Screen_Run.get(RunID)
    valLog = Validation_Log("Gen_Masterplates")

    logNumbers = {'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                  'Valid Plates':0, 'Rejected Plates':0, 'Failed Plates':0,
                  'Processed AssayData':0,'Uploaded AssayData':0,
                  'Inhibition AssayData':0, 'MIC AssayData':0, 'CC50 AssayData':0, 'HC50 AssayData':0, 
                  'Empty':0}

    if nRacks > 0:
        Barcodes = get_BarcodeScans(DirName, RackFileList, valLog=valLog)
    
    if nFiles > 0 :
        dfMP = gen_Motherplates_PSPrep(DirName,PrepFileList[0],Barcodes,valLog=valLog,verbose=0)

        if generate and len(dfMP)>0:
            xlsFile = f"{RunID}_PS_Motherplates.xlsx"
            print(f" Generate Download {DirName} {xlsFile}")
            dfMP.to_excel(os.path.join(DirName,xlsFile))

            # req = HttpResponse(content_type='application/vnd.ms-excel')
            # req['Content-Disposition'] = f'attachment; filename={_xls_name}'
            # dfMP.to_excel(req)


    valLog.select_unique()
    
    return(valLog)
