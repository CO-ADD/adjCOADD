import os
import pandas as pd
from decimal import Decimal

from applib.process.process_progress import Process_Progress
from dsample.models import Project
from dcollab.models import Collab_Group, Collab_User, Collab_Membership, Organisation
from dplate.models import MasterPlate, MasterWell


from applib.plate.stockprep import  read_Stock_Prepsheet_XLS
from applib.project.import_project import get_CompoundSubmisssion_xlsx, parse_SampleInfo_Sheet, parse_ContactInfo_Sheet, Upload_Project_Collab
from dsample.utils.summary import update_project_summary
from dsample.utils.compounds import Upload_COADD_Compound

from django.conf import settings
import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------
class Load_Project_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    process_name = 'Upload_Project'

    def __init__(self, request, dirname, file_list, project_id=None,
                 upload_content={'Samples':True,'Contacts':True},
                 upload=False, overwrite=False, appuser=None,
                 progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.project_id = project_id
        self.upload_content = upload_content
        self.upload = upload
        self.overwrite = overwrite

    def run(self):
        nFiles = self.n_files
        nUploads = 0

        valLog = self.valLog
        procDict = {}

        if nFiles > 0:
            for i in range(nFiles):
                self.report_progress(i + 1, nFiles, self.file_list[i], 'Reading Compound Submission')

                if settings.DEBUG:
                    print(f" [Upload_Project] {i+1:3d}/{nFiles:3d} - {self.file_list[i]}  [{self.appuser}] ")

                _dictSheets = get_CompoundSubmisssion_xlsx(os.path.join(self.dirname,self.file_list[i]), valLog=valLog)

                # - Project Init ------------------------------------------------------
                djPrj = None
                if self.project_id:
                    djPrj = Project.get(self.project_id)
                    if djPrj is None:
                        djPrj = Project()
                        djPrj.project_id = self.project_id
                        valLog.add_warning(f"New Project",f"{self.project_id}",f"{self.project_id} does not exist")
                        self.upload_content['Contacts'] = True
                    else:
                        valLog.add_info(f"Existing Project",f"{djPrj}")
                else:
                    djPrj = Project()
                    valLog.add_warning(f"New Project","<New ID>","ProjectID not provided")
                    self.upload_content['Contacts'] = True

                if _dictSheets['Contacts'] is not None:
                    _Contacts,_PrjTitle = parse_ContactInfo_Sheet(_dictSheets['Contacts'],valLog=valLog)
                    djPrj.project_name = _PrjTitle

                if _dictSheets['Samples'] is not None:
                    _Samples = parse_SampleInfo_Sheet(_dictSheets['Samples'],valLog=valLog)

                # - Contacts and Project ----------------------------------------------
                if self.upload_content['Contacts'] and _dictSheets['Contacts'] is not None:
                    djPrj.project_name = _PrjTitle
                    Upload_Project_Collab(djPrj, _Contacts, upload=self.upload, overwrite=self.overwrite, valLog=valLog)

                #- Samples ----------------------------------------------
                if self.upload_content['Samples'] and _dictSheets['Samples'] is not None:
                    _nSamples = len(_Samples)
                    _done = 0
                    for _smp in _Samples:
                        _done += 1
                        _smp['project_id'] = str(djPrj)
                        self.report_progress(_done, _nSamples, _smp.get('compound_id', ''), 'Uploading sample')
                        Upload_COADD_Compound(_smp, upload=self.upload, overwrite=self.overwrite, valLog=valLog)

                procDict['pk'] = djPrj.project_id

        else:
            print(f" [Upload_Project] No Xlsx to process in {self.dirname}  ")

        self.finish()
        valLog.select_unique()
        print(f" [Load_Project_Process] : {procDict} ")
        return(valLog, procDict)

#-----------------------------------------------------------------------------------
class Summary_Project_Process(Process_Progress):
#-----------------------------------------------------------------------------------
    process_name = 'Summary_Project'

    def __init__(self, request, project_id, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, appuser=appuser, progress_cb=progress_cb, progress_key=progress_key)
        self.project_id = project_id

    def run(self):
        self.report_progress(0, 1, self.project_id, 'Updating Project summary ...')
        djPrj = Project.get(self.project_id)
        print(f" [Project] Update Summary [{djPrj}]")
        update_project_summary(djPrj)
        djPrj.save()
        self.finish()
        return self.valLog

#-----------------------------------------------------------------------------------
class Upload_StockPrep_Process(Process_Progress):
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
    process_name = 'Upload_StockPrep'

    def __init__(self, request, dirname, file_list, project_id=None, upload=False,
                 overwrite=False, appuser=None, progress_cb=None, progress_key=None):
        super().__init__(request, dirname, file_list, appuser=appuser,
                         progress_cb=progress_cb, progress_key=progress_key)
        self.project_id = project_id
        self.upload = upload
        self.overwrite = overwrite

    def run(self):
        nFiles = self.n_files
        nUploads = 0

        djPrj = Project.get(self.project_id)
        valLog = self.valLog

        if nFiles > 0:
            for i in range(nFiles):
                valLog.add_info('Read StockPrep File', self.file_list[i])
                self.report_progress(i + 1, nFiles, self.file_list[i], 'Reading StockPrep file')

                if settings.DEBUG:
                    print(f" [Upload_StockPrep] {i+1:3d}/{nFiles:3d} - {self.file_list[i]}  [{djPrj}]  [{self.appuser}] ")
                lstMP = read_Stock_Prepsheet_XLS(os.path.join(self.dirname,self.file_list[i]),valLog=valLog)

                print(lstMP)
                nPlates = len(lstMP)
                _done_plates = 0
                # Process MP's
                for _mpid in lstMP:
                    _done_plates += 1
                    self.report_progress(_done_plates, nPlates, str(lstMP[_mpid]['plate']), 'Validating motherplate')
                    print(f" [Upload_StockPrep] Validating: {lstMP[_mpid]['plate']} ")
                    validStatus = True
                    validDict = {}

                    lstMP[_mpid]['plate'].set_defaults_model()

                    validDict = lstMP[_mpid]['plate'].validate_model(WellData=True, verbose = 0)
                    if validDict:
                        validStatus = False
                        for c in validDict:
                            print(f" [Upload_StockPrep] validDict: {c} ")

                    if self.upload and validStatus:
                        if lstMP[_mpid]['new'] or lstMP[_mpid]['empty'] or self.overwrite:
                            print(f" [Upload_StockPrep] Saving: {lstMP[_mpid]['plate']} [Overwrite: {self.overwrite}]")
                            lstMP[_mpid]['plate'].save(verbose=0)
                            nUploads += 1
            if self.upload:
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
            print(f"[Upload_StockPrep] No Xlsx to process in {self.dirname}  ")

        self.finish()
        valLog.select_unique()
        return valLog




