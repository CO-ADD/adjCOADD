import os

from applib.logging.validation_log import Validation_Log
from applib.process.process_progress import Process_Progress

from dscreen.models import Screen_Run
from dplate.models import Labware, TestPlate, TestWell, MasterPlate
from dorganism.models import Organism_Batch
from dgene.models import Genome_Sequence

from applib.plate.multimode_reader import multimodereader_xls
#from applib.plate.masterplates import read_motherplate_prepsheet_xls
from applib.plate.plateprep import read_Motherplates_Prepsheet_XLS, read_TestPlateList_Prepsheet_XLS, get_BarcodeScans, gen_Motherplates_PSPrep, gen_Motherplates_HCPrep
from applib.plate.testplates import add_mother_to_testplate
from applib.bio.doseresponse import process_testplate_doseresponse
from applib.data.dfutils import get_Xlxs_Sheets
from applib.data.set_fielddata import set_model_from_dict
from dscreen.utils.summary import update_screenrun_summary
from adjcoadd.constants import DR_CLASSES

import pandas as pd

from django.conf import settings
import logging
logger = logging.getLogger(__name__)


#-----------------------------------------------------------------------------------
class Summary_ScreenRun_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    process_name = 'Summary_ScreenRun'

    def __init__(self, request, run_id, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, appuser=appuser, progress_cb=progress_cb, progress_key=progress_key)
        self.run_id = run_id

    def run(self):
        self.report_progress(0, 1, self.run_id, 'Updating ScreenRun summary ...')
        djRun = Screen_Run.get(self.run_id)
        print(f" [ScreenRun] Update Summary [{djRun}]")
        update_screenrun_summary(djRun)
        djRun.save()
        self.finish()
        return self.valLog


#-----------------------------------------------------------------------------------
class Upload_ReadOuts_Process(Process_Progress):
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
    process_name = 'Upload_ReadOuts'

    def __init__(self, request, dirname, file_list, run_id=None, upload=False,
                 overwrite=False, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.run_id = run_id
        self.upload = upload
        self.overwrite = overwrite

    def run(self):
        nFiles = self.n_files
        nUploads = 0

        djRun = Screen_Run.get(self.run_id)
        valLog = self.valLog

        if nFiles > 0:
            for i in range(nFiles):
                valLog.add_info('Read Readout File', self.file_list[i])
                self.report_progress(0, 1, self.file_list[i], f'Reading file {i+1}/{nFiles} ...')

                if settings.DEBUG:
                    print(f" [Upload_ReadOuts] {i+1:3d}/{nFiles:3d} - {self.file_list[i]}  [{djRun}]  [{self.appuser}] ")
                lstTP = multimodereader_xls(os.path.join(self.dirname,self.file_list[i]),valLog=valLog)

                nPlates = len(lstTP)
                _progress_plates = 0
                
                for _tp in lstTP:
                    _progress_plates += 1
                    self.report_progress(_progress_plates, nPlates, _tp['plate_id'], 'Reading TestPlate')
                    
                    if settings.DEBUG:
                        logger.info(f" [Upload_ReadOuts] Validating: {_tp['plate']} {_progress_plates}/{nPlates}")

                    validStatus = True
                    validDict = {}
                    _tp['plate'].run_id = djRun
                    _tp['plate'].set_defaults_model()
                    validDict = _tp['plate'].validate_model(WellData=False, verbose = 0)
                    if validDict:
                        validStatus = False
                        for c in validDict:
                            logger.error(f" [Upload_ReadOuts] validDict: {c} ")

                    if self.upload and validStatus:
                        if _tp['new'] or self.overwrite:
                            _tp['plate'].save(verbose=0)
                            logger.info(f" [Upload_ReadOuts] Saving: {_tp['plate_id']} [Overwrite: {self.overwrite}]")
                            nUploads += 1

            if self.upload:
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
            valLog.add_error(f"No Xlsx to process",f" in {self.dirname} ")

        self.finish()
        valLog.select_unique()
        return valLog

#-----------------------------------------------------------------------------------
class Upload_Motherplates_Process(Process_Progress):
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
    process_name = 'Upload_Motherplates'

    def __init__(self, request, dirname, file_list, run_id=None, upload=False,
                 overwrite=False, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.run_id = run_id
        self.upload = upload
        self.overwrite = overwrite

    def run(self):
        nFiles = self.n_files
        nUploads = 0

        djRun = Screen_Run.get(self.run_id)
        valLog = self.valLog

        if nFiles > 0:
            for i in range(nFiles):
                valLog.add_info('Read PlatePrep File', self.file_list[i],"[MotherPlates]")
                self.report_progress(i + 1, nFiles, self.file_list[i], 'Reading PlatePrep file')

                print(f" [Upload_Motherplates] {i+1:3d}/{nFiles:3d} - {self.file_list[i]}  [{djRun}]  [{self.appuser}] ")
                lstMP = read_Motherplates_Prepsheet_XLS(os.path.join(self.dirname,self.file_list[i]),valLog=valLog)

                nPlates = len(lstMP)
                _done_plates = 0
                for _mp in lstMP:
                    _done_plates += 1
                    self.report_progress(_done_plates, nPlates, str(_mp['plate']), 'Validating motherplate')
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

                    if self.upload and validStatus:
                        if _mp['new'] or self.overwrite:
                            _mp['plate'].save(verbose=0)
                            print(f" [Upload_MotherPlates] Saving: {_mp['plate']} [Overwrite: {self.overwrite}]")
                            nUploads += 1
            if self.upload:
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
            print(f"[Upload_Motherplates] No Xlsx to process in {self.dirname}  ")

        self.finish()
        valLog.select_unique()
        return valLog
#-----------------------------------------------------------------------------------
class Upload_TestplateList_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    """
    Progress-reporting version of the TestplateList upload.

    Usage (explicit)::

        process = Upload_TestplateList_Process(
            request, dirname, file_list, run_id=pk,
            upload=True, apply_mp=True, appuser=request.user,
        )
        valLog = process.run()

    or, function-style::

        valLog = Upload_TestplateList_Process.execute(
            request, dirname, file_list, run_id=pk,
            upload=True, apply_mp=True, appuser=request.user,
        )

    ``progress_cb`` receives dicts of the form
    ``{'done', 'total', 'percent', 'item', 'message'}``; when omitted, progress
    is only written to the cache under ``progress_key`` for polling.
    """

    process_name = 'Upload_TestplateList'

    def __init__(self, request, dirname, file_list, run_id=None,
                 upload=False, apply_mp=False, overwrite=False, appuser=None,
                 progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list,
                         appuser=appuser, progress_cb=progress_cb, progress_key=progress_key)

        self.run_id = Screen_Run.get(run_id)
        self.upload = upload
        self.apply_mp = apply_mp
        self.overwrite = overwrite
        
        self.tp_lst = {}
        self.mp_lst = {}
        self.ass_lst = {}

    # ----------------------------------------------------------------
    def _process_files(self):
    # ----------------------------------------------------------------
        for i in range(self.n_files):
            self.valLog.add_info('Read PlatePrep File', self.file_list[i],"[TestPlateList, Assay]")
            
            self.report_progress(i, self.n_files, message= f'Reading file {i+1}/{self.n_files} ...')
            if settings.DEBUG:
                print(f" [Upload_TestplateList] {i+1:3d}/{self.n_files:3d} - {self.file_list[i]}  [{self.run_id}]  [{self.appuser}] ")

            dictTP,dictMP,dictAss = read_TestPlateList_Prepsheet_XLS(os.path.join(self.dirname,self.file_list[i]),
                                                            apply_mp=self.apply_mp,valLog=self.valLog)
            self.tp_lst.update(dictTP)
            self.mp_lst.update(dictMP)
            self.ass_lst.update(dictAss)
            
            for mp in self.mp_lst:
                self.mp_lst[mp] = MasterPlate.get(mp,WellData=True)
                if self.mp_lst[mp] is None:
                    self.valLog.add_error('MotherPlate not found ', mp, "Add MNotherPlate or check ID")
                else:
                    if str(self.mp_lst[mp].run_id) != str(self.run_id):
                        self.valLog.add_warning("Different RunID",f" {str(self.mp_lst[mp].plate_id)}",f"From {str(self.mp_lst[mp].run_id)}","")                     

        self.valLog.add_info('PlatePrep ', f"Testplates  : {len(self.tp_lst)}","") 
        self.valLog.add_info('PlatePrep ', f"Motherplates: {len(self.mp_lst)}","") 
        self.valLog.add_info('PlatePrep ', f"Assays      : {len(self.ass_lst)}","") 

    # ----------------------------------------------------------------
    def _process_testplates(self):
    # ----------------------------------------------------------------
        progressCount = 0
        progressTotal = len(self.tp_lst)
        
        progressLog = {'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                      'Valid Plates':0, 'Rejected Plates':0, 'Failed Plates':0,
                      'Processed AssayData':0,'Uploaded AssayData':0,
                      'Inhibition AssayData':0, 'MIC AssayData':0, 'CC50 AssayData':0, 'HC50 AssayData':0,
                      'Empty':0}
    
        for key in self.tp_lst:
            progressCount += 1
            _tp = self.tp_lst[key]['plate']

            self.report_progress(progressCount, progressTotal, _tp.plate_id, f'Processing {_tp.plate_id}')
            
            if settings.DEBUG:
                print(f" [Upload_TestplateList] {progressCount:3d}/{progressTotal:3d} - {key}  [{self.run_id}]  ")
            
            validDict = _tp.validate_model(WellData=False, verbose = 0)
            if validDict:
                validStatus = False
                for c in validDict:
                    print(f" [Upload_TestPlateList] validDict: {c} ")

            if self.apply_mp:
                
                # Apply MotherPlates -> Add Compounds
                # -------------------------
                if hasattr(_tp,'motherplate_ids'):
                    if settings.DEBUG:
                        print(f" [Upload_TestplateList] {_tp.plate_id} Apply MP: {_tp.motherplate_ids}")

                    _tp.clear_cmpbatch_data()
                    mp_ids = getattr(_tp,'motherplate_ids')
                    for mp in mp_ids:
                        if mp != '-':
                            print(f"add {mp} {self.mp_lst[mp]} -> {_tp}")
                            add_mother_to_testplate(self.mp_lst[mp],_tp)
                            print('done')
                    _tp.update_n('n_samples')
                    
                else:
                    validStatus = False
                    self.valLog.add_error("No MotherPlates",key,"No MotherPlate_IDs","Correct TestPlateList")

                print(f"{validStatus} and {_tp.n_wells} and {_tp.n_reads} and {_tp.control_layout}")
                # Analyze Testplate
                # -----------------
                if validStatus and _tp.n_wells > 0 and _tp.n_reads > 0 and _tp.control_layout:

                    if settings.DEBUG:
                        print(f" [Upload_TestplateList] {_tp.plate_id} Apply Layout: {_tp.control_layout}")
                    _n_layout = _tp.apply_layout()
                    
                    if _n_layout > 0:
                        # Analyze Testplate -> Inhibition
                        # -------------------------------
                        if settings.DEBUG:
                            print(f" [Upload_TestPlateList] {_tp.plate_id} Apply Calc ")

                        _tp.calc_inhibition()
                        
                        if settings.DEBUG:
                            print(f" [Upload_TestPlateList] {_tp.plate_id} {_tp.plate_quality} {_tp.zfactor} {_tp.n_inhibitions}")

                        _dr_list = []
                        if str(_tp.plate_quality) == 'Valid' and str(_tp.result_type) in DR_CLASSES:
                            if settings.DEBUG:
                                print(f" [Upload_TestPlateList] {_tp.plate_id} ApplyDoseReponse  ")

                            _dr_list = process_testplate_doseresponse(_tp)

                            if settings.DEBUG:
                                print(f" [Upload_TestplateList] {key} DR: {_tp.result_type} {len(_dr_list)}")

                        elif str(_tp.plate_quality) == 'Rejected':
                            self.valLog.add_error("Testplate rejected",_tp.plate_id,f"Zfactor : {_tp.zfactor}")
                else:
                    validStatus = False
                    self.valLog.add_error("Calc Error",_tp.plate_id,"Either no Wells, Readout or Layout","Correct TestPlateList")

            if self.upload and validStatus:
                if settings.DEBUG:
                    print(f" [Upload_TestPlateList] {key} Updating: {_tp.plate_quality} {_tp.zfactor} {_tp.n_inhibitions}")
                _tp.save()
                
                for _dr in _dr_list:
                    _dr.save_assaydata(overwrite=True)
                            
    # ----------------------------------------------------------------
    def run(self):
        self._process_files()
        self._process_testplates()
        
        self.valLog.select_unique()
        return self.valLog
        
    # ----------------------------------------------------------------
    def run_old(self):
        nFiles = self.n_files
        nUploads = 0

        djRun = Screen_Run.get(self.run_id)
        valLog = self.valLog

        progressCount = 0
        progressTotal = 0
        progressLog = {'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                      'Valid Plates':0, 'Rejected Plates':0, 'Failed Plates':0,
                      'Processed AssayData':0,'Uploaded AssayData':0,
                      'Inhibition AssayData':0, 'MIC AssayData':0, 'CC50 AssayData':0, 'HC50 AssayData':0,
                      'Empty':0}
        
        if nFiles > 0:
            _MPS = {}
            for i in range(nFiles):
                valLog.add_info('Read PlatePrep File', self.file_list[i],"[TestPlateList, Assay]")
                self.report_progress(0, 1, self.file_list[i], f'Reading file {i+1}/{nFiles} ...')
                if settings.DEBUG:
                    print(f" [Upload_TestplateList] {i+1:3d}/{nFiles:3d} - {self.file_list[i]}  [{djRun}]  [{self.appuser}] ")

                lstTP,lstAss = read_TestPlateList_Prepsheet_XLS(os.path.join(self.dirname,self.file_list[i]),apply_mp=self.apply_mp,valLog=valLog)

                if valLog.if_noerrors():
                    
                    nPlates = len(lstTP)
                    progressTotal = nPlates
                    progressCount = 0
                    
                    for _tpid in lstTP:
                        progressCount += 1
                        _tp = lstTP[_tpid]['plate']
                        
                        validStatus = True
                        self.report_progress(progressCount, progressTotal, _tp.plate_id, f'Processing {_tp.plate_id}')
                        if settings.DEBUG:
                            print(f" [Upload_TestplateList] {_tp.plate_id} Validating:  {_tp.result_type} {_tp.assay_id} [{progressCount}/{progressTotal}]")
                            
                        validDict = _tp.validate_model(WellData=False, verbose = 0)
                        if validDict:
                            validStatus = False
                            for c in validDict:
                                print(f" [Upload_TestPlateList] validDict: {c} ")

                        # Get and Check Motherplates
                        if hasattr(_tp,'motherplate_ids'):
                            mp_ids = getattr(_tp,'motherplate_ids')
                            for mp in mp_ids:
                                if mp != '-':
                                    if mp not in _MPS:
                                        _MPS[mp] = MasterPlate.get(mp,WellData=True)
                                        if str(_MPS[mp].run_id) != self.run_id:
                                            valLog.add_warning("Different RunID",f" {str(_MPS[mp].plate_id)}",f"From {str(_MPS[mp].run_id)}","")                            
                        
                        if self.apply_mp:
                            # Apply MP -> Add Compounds
                            # -------------------------
                            if hasattr(_tp,'motherplate_ids'):
                                if settings.DEBUG:
                                    print(f" [Upload_TestplateList] {_tp.plate_id} Apply MP: {_tp.motherplate_ids}")

                                _tp.clear_cmpbatch_data()
                                mp_ids = getattr(_tp,'motherplate_ids')
                                for mp in mp_ids:
                                    if mp != '-':
                                        # if mp not in _MPS:
                                        #     _MPS[mp] = MasterPlate.get(mp,WellData=True)
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
                                        progressLog['Valid Plates'] += 1

                                        if str(_tp.result_type) in DR_CLASSES:
                                            # Analyze Testplate -> DR AssayData
                                            # ---------------------------------
                                            _dr_list = process_testplate_doseresponse(_tp)

                                            # Analyze Testplate -> AssayData
                                            # ------------------------------
                                            progressLog['Processed AssayData'] += len(_dr_list)
                                            for _dr in _dr_list:
                                                progressLog[f'{_dr.dr_type} AssayData'] += 1

                                            if settings.DEBUG:
                                                print(f" [Upload_TestplateList] {_tpid} DR: {_tp.result_type} {len(_dr_list)}")

                                        elif str(_tp.result_type) == 'Inhibition':
                                            progressLog[f'Inhibition AssayData'] += _tp.n_inhibitions

                                            if settings.DEBUG:
                                                print(f" [Upload_TestplateList] {_tpid} SC: {_tp.result_type} {_tp.n_inhibitions}")
                                        else:
                                            if settings.DEBUG:
                                                print(f" [Upload_TestplateList] {_tpid} ??: {_tp.result_type}")

                                    elif str(_tp.plate_quality) == 'Rejected':
                                        progressLog['Rejected Plates'] += 1
                                    else:
                                        progressLog['Failed Plates'] += 1

                            else:
                                validStatus = False
                                valLog.add_error("Calc Error",_tp.plate_id,"Either no Wells, Readout or Layout","Correct TestPlateList")

                        if self.upload and validStatus:
                            if settings.DEBUG:
                                print(f" [Upload_TestPlateList] {_tpid} Updating: {_tp.plate_quality} {_tp.zfactor} {_tp.n_inhibitions}")
                            _tp.save()
                            for _dr in _dr_list:
                                _dr.save_assaydata(overwrite=True)
                            progressLog['Uploaded Plates'] += 1
                            progressLog['Uploaded AssayData'] += len(_dr_list)

                    # Message valLog
                    if self.apply_mp:

                        #TestPlates
                        valLog.add_info('TestPlates',f"Valid: {progressLog['Valid Plates']} ","Testplates pass QC","Confirm Testplate")
                        if progressLog['Rejected Plates'] > 0:
                            valLog.add_warning('TestPlates',f"Rejected: {progressLog['Rejected Plates']} ","Testplates Rejected by ZFactor QC","Check Testplate/Layout")
                        else:
                            valLog.add_info('TestPlates',f"Rejected: {progressLog['Rejected Plates']} ","Testplates Rejected by ZFactor QC")

                        if progressLog['Failed Plates'] > 0:
                            valLog.add_warning('TestPlates',f"Failed: {progressLog['Failed Plates']} ","Testplates fail ","Check Testplate")
                        else:
                            valLog.add_info('TestPlates',f"Failed: {progressLog['Failed Plates']} ","Testplates fail ")

                        # AssayData
                        valLog.add_info('SC AssayData',f"Inhibition: {progressLog['Inhibition AssayData']} ","Inhibition AssayData generated ")
                        valLog.add_info('DR AssayData',f"Processed: {progressLog['Processed AssayData']} ","Dosereponse AssayData generated ")
                        valLog.add_info('DR AssayData',f"MIC: {progressLog['MIC AssayData']} ","MIC AssayData generated ")
                        valLog.add_info('DR AssayData',f"CC50: {progressLog['CC50 AssayData']} ","CC50 AssayData generated ")
                        valLog.add_info('DR AssayData',f"HC50: {progressLog['HC50 AssayData']} ","HC50 AssayData generated ")
                        if progressLog['Processed AssayData'] > 0:
                            valLog.add_warning('TestPlates',f"Failed: {progressLog['Failed Plates']} ","Testplates fail ","Check Testplate")
                        else:
                            valLog.add_info('TestPlates',f"Failed: {progressLog['Failed Plates']} ","Testplates fail ")

                    else:
                        valLog.add_warning('TestPlates',f"No Analysis",f"No Compound, Layout or Analysis applied","Select Apply_MP")

                    if self.upload:
                        valLog.add_info('TestPlates',f"Uploaded: {progressLog['Uploaded Plates']} ","Testplates uploaded")
                        _not_uploaded = progressLog['Processed Plates']- progressLog['Uploaded Plates']
                        if  _not_uploaded > 0:
                            valLog.add_warning('TestPlates',f"Not Uploaded: {_not_uploaded} ","Testplates not uploaded","Check Testplate")
                        valLog.add_info('DR AssayData',f"Uploaded: {progressLog['Uploaded AssayData']} ","Dosereponse AssayData uploaded")

                    else:
                        valLog.add_warning('TestPlates',f"Uploaded: 0",f"No Testplate uploaded to database","Select Upload")

        else:
            print(f"[Upload_TestplateList] No Xlsx to process in {self.dirname}  ")

        # Signal completion
        self.finish()
        valLog.select_unique()
        return valLog


#-----------------------------------------------------------------------------------
class Gen_Masterplates_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    """
    Generates Motherplates (PSR/HCR) from a PrepSheet and TubeRack scans.

        request        : Request object (user, etc.); may be None
        dirname        : Folder containing the files
        prep_file_list : [PrepSheet.xlsx] - the PrepSheet workbook
        rack_file_list : [TubeRack.csv, ...] - TubeRack scans
        run_id         : Screen_Run ID
        generate       : Write the generated Motherplate xlsx when True
        appuser        : User instance

    Progress is reported at task level:
      (1) scanning TubeRacks, (2) generating the Motherplate.
    """

    process_name = 'Gen_Masterplates'

    def __init__(self, request, dirname, prep_file_list, rack_file_list, run_id=None,
                 generate=False, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, dirname, prep_file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.prep_file_list = prep_file_list
        self.rack_file_list = rack_file_list
        self.run_id = run_id
        self.generate = generate

    def run(self):
        nFiles = len(self.prep_file_list) if self.prep_file_list else 0
        nRacks = len(self.rack_file_list) if self.rack_file_list else 0

        djRun = Screen_Run.get(self.run_id)
        valLog = self.valLog

        Barcodes = {}
        if nRacks > 0:
            self.report_progress(0, 2, 'Scanning TubeRacks', f'Scanning {nRacks} rack file(s)')
            Barcodes = get_BarcodeScans(self.dirname, self.rack_file_list, valLog=valLog)

        if djRun is not None and nFiles > 0:
            dfMP = pd.DataFrame()
            self.report_progress(1, 2, 'Generating Motherplates', 'Generating Motherplates from PrepSheet')

            if str(djRun.run_type) == 'PSR':
                dfMP = gen_Motherplates_PSPrep(self.dirname, self.prep_file_list[0], Barcodes, valLog=valLog, verbose=0)
                xlsFile = f"{self.run_id}_PSPrep_Motherplates.xlsx"

            elif str(djRun.run_type) == 'HCR':
                dfMP = gen_Motherplates_HCPrep(self.dirname, self.prep_file_list[0], Barcodes, valLog=valLog, verbose=0)
                xlsFile = f"{self.run_id}_HCPrep_Motherplates.xlsx"

            if self.generate and len(dfMP) > 0:
                dfMP.to_excel(os.path.join(self.dirname, xlsFile))

        self.finish()
        valLog.select_unique()
        return valLog

#-----------------------------------------------------------------------------------
class Upload_Sequences_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) genome sequences from a Sequencing Excel file.

        request    : Request object (user, etc.); may be None
        dirname    : Folder containing the file
        file_list  : [Sequencing.xlsx] - the Sequencing workbook
        run_id     : Screen_Run ID
        upload     : Validation only (False) or Validation and Upload (True)
        overwrite  : On Upload overwrite existing data
        appuser    : User instance

    Progress is reported per sequence row processed.
    """

    process_name = 'Upload_Sequences'

    def __init__(self, request, dirname, file_list, run_id=None, sheets=None,
                 upload=False, overwrite=False, appuser=None,
                 progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.run_id = run_id
        self.sheets = sheets
        self.upload = upload
        self.overwrite = overwrite

    def run(self):
        # Xlsx Definitions
        SeqRun_Sheets = {'Seq': None,}
        SeqRun_List_Fields = ['runsample_file','runsample_dir','runsample_name',]
        SeqRun_List_Dicts = ['seq_type','seq_method']
        unqSeqCodes = set()

        nFiles = self.n_files

        djRun = Screen_Run.get(self.run_id)
        valLog = self.valLog

        if nFiles > 0:
            get_Xlxs_Sheets(self.dirname, self.file_list[0], SeqRun_Sheets, FillNA=None, valLog=valLog)

            df = SeqRun_Sheets['Seq']
            df = df[df['run_id'] == str(djRun)]
            if len(df) > 0:
                nSeq = 0
                nRows = len(df)
                done = 0
                for idx, row in df.iterrows():
                    done += 1
                    _item = row.get('seq_code')
                    if pd.isna(_item):
                        _item = row.get('orgbatch_id', '')
                    self.report_progress(done, nRows, str(_item), 'Processing sequence')

                    # Check OrgBatch
                    djOrgBatch = Organism_Batch.get(row['orgbatch_id'])
                    if djOrgBatch is None:
                        valLog.add_error('OrgBatch not found', row['orgbatch_id'], f"Correct OrgBatch_ID or Import new OrgBatch")
                    else:
                        nSeq += 1

                    # SeqCode - Check Duplicates
                    if pd.isna(row['seq_code']):
                        row['seq_code'] = Genome_Sequence.gen_seq_code(row['orgbatch_id'], str(djRun))

                    if row['seq_code'] in unqSeqCodes:
                        valLog.add_error('Duplicate SeqCode', f"{row['seq_code']} - {row['orgbatch_id']} {str(djRun)} ", f"Correct SeqCode with unique values")
                    else:
                        unqSeqCodes.add(row['seq_code'])

                    # New Seq - Check SeqCode
                    djSeq = Genome_Sequence.get(SeqCode=row['seq_code'])
                    if djSeq is None:
                        djSeq = Genome_Sequence()
                        djSeq.orgbatch_id = djOrgBatch
                        djSeq.run_id = djRun
                        djSeq.seq_code = row['seq_code']
                    else:
                        valLog.add_error('SeqCode exists already', f"{row['seq_code']} ", f"Correct SeqCode with unique values")

                    # Set SeqAttributes
                    validStatus = set_model_from_dict(djSeq, row,
                                                      list_Fields=SeqRun_List_Fields,
                                                      list_Dicts=SeqRun_List_Dicts,
                                                      valLog=valLog)
                    # Uploading
                    if self.upload and validStatus:
                        if settings.DEBUG:
                            print(f" [Upload Sequence] {djSeq} ")
                        djSeq.save()

                valLog.add_info('Reading SeqRun file', f"{nSeq} Sequences", "")
            else:
                valLog.add_error('[Seq] no sequences for RunID', str(self.run_id), f"Correct Run_ID in [Seq].SeqRun_ID")

        self.finish()
        return valLog

