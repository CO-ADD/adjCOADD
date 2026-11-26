import os, sys
import pandas as pd
import logging
from tqdm import tqdm
from decimal import Decimal

from adjcoadd import settings
from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dplate.models import TestPlate, MasterPlate
from dscreen.models import Screen_Run
from adjcoadd.constants import DR_CLASSES, SYN_CLASSES
from applib.plate.testplates import add_mother_to_testplate
from applib.bio.doseresponse import calc_testplate_doseresponse, reset_testplate_assaydata
from applib.bio.synergy import calc_testplate_synergy


logger = logging.getLogger(__name__)

#================================================================================================
class Command(BaseCommand):
#================================================================================================
    
    command_process = "Process TestPlates"

    def add_arguments(self, parser):
        # parser.add_argument("--projectid")
        # parser.add_argument("--xlsx_file")
        parser.add_argument("--runid",action="store",default=None,help='RunID')
        parser.add_argument("--plateid",action="store",default=None,help='PlateID')
        
        parser.add_argument("--apply",action="store_true",default=False)
        
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        parser.add_argument("--debug",action="store_true",default=False)
        
    
    #------------------------------------------------------------------  
    def _process_testplate(self, PlateID, *args, **kwargs):
    #------------------------------------------------------------------  
        
        apply = kwargs.get('apply', False)
        upload = kwargs.get('upload', False)
        overwrite = kwargs.get('overwrite', False)
        
        verbose = kwargs.get('verbose', 0)
        valLog = kwargs.get('valLog', None)
        
        _Valid = True
        djTP = TestPlate.get(PlateID)
        
        if djTP is None:
            logger.error(f" [{self.command_process}] TestPlateID {PlateID} not found")
            return()
        
        if apply:
            # Apply MotherPlates
            #-----------------------------------------------
            if hasattr(djTP,'motherplate_ids'):
                if verbose > 0:
                    logger.debug(f" [{self.command_process}] Testplate Apply MotherPlates")
                djTP.clear_cmpbatch_data()
                mp_ids = getattr(djTP,'motherplate_ids')
                for mp in mp_ids:
                    if mp != '-':
                        djMP = MasterPlate.get(mp)
                        if djMP is None:
                            _Valid = False
                            print(f" [{self.command_process}] MasterPlateID {mp} not found")
                            if valLog:
                                valLog.add_error("No MotherPlates",mp,"No MotherPlate_IDs","Correct TestPlateList")
                        else:                            
                            add_mother_to_testplate(djMP,djTP)      
                djTP.update_n('n_samples')
            
            # Apply Layout
            #-----------------------------------------------
            if _Valid and djTP.n_wells > 0 and djTP.n_reads > 0 and djTP.control_layout:
                if verbose > 0:
                    logger.debug(f" [{self.command_process}] Testplate Apply Layout")
                djTP.n_layout = djTP.apply_layout()            
                
            # Calc Inhibition
            #-----------------------------------------------
            if _Valid and djTP.n_layout > 0:
                if verbose > 0:
                    logger.debug(f" [{self.command_process}] Testplate Calculate Inhibition")
                djTP.calc_inhibition()
                
                if str(djTP.plate_quality) == 'Rejected':
                    if valLog:
                        valLog.add_warning("Testplate rejected",djTP.plate_id,f"Zfactor : {djTP.zfactor}")

        
        if _Valid and upload and overwrite:
            # Reset AssayData - MIC, CC50, HC50
            #-----------------------------------------------
            if str(djTP.result_type) in DR_CLASSES:
                if verbose > 0:
                    logger.debug(f" [{self.command_process}] Reset AssayData_{str(djTP.result_type)}")
                reset_testplate_assaydata(djTP, upload=upload, overwrite=overwrite, verbose=verbose)
                
        if _Valid:
            # Process TestPlate -> AssayData
            #-----------------------------------------------
            if str(djTP.result_type) in DR_CLASSES:
                _dr_list = calc_testplate_doseresponse(djTP,upload, overwrite=overwrite,verbose=verbose)
            elif str(djTP.result_type) in SYN_CLASSES:
                _dr_dict, _iso_lst, _syn_lst = calc_testplate_synergy(djTP,upload=upload, overwrite=overwrite, verbose=verbose)
            
                
    #------------------------------------------------------------------  
    def handle(self, *args, **options):
    #------------------------------------------------------------------  

        if options['debug']:
            logger.setLevel(logging.DEBUG)
            verbose=2
        else:
            logger.setLevel(logging.INFO)
            verbose=0
            
        #valLog = Validation_Log("Process_TestPlates")
        
        _dr_list = []
        if options['plateid']:
            logger.info(f" [{self.command_process}] Processing TestPlateID {options['plateid']} ...")
            self._process_testplate(options['plateid'], upload=options['upload'],overwrite=options['overwrite'],apply=options['apply'],verbose=verbose)
                                        
        elif options['runid']:
            djRun = Screen_Run.get(options['runid'])
            if djRun:
                qryTP = TestPlate.objects.filter(run_id=djRun).values_list('plate_id')
                nTP = qryTP.count()
                logger.info(f" [{self.command_process}] Processing RunID {djRun} with {nTP} testplates ...")

                for _tp_id in tqdm(qryTP, total=nTP, desc='TestPlates'):
                    self._process_testplate(_tp_id[0], upload=options['upload'],overwrite=options['overwrite'],apply=options['apply'],verbose=verbose)
            else:
                logger.error(f" [{self.command_process}] RunID {options['runid']} NOT Found")
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

