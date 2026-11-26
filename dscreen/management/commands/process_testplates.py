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
    def _process_testplate(self, djTP, *args, **options):
    #------------------------------------------------------------------  
        
        _Valid = True
        djTP = TestPlate.get(options['plateid'])
        
        if djTP is None:
            print(f" [{self.command_process}] TestPlateID {options['plateid']} not found")
            return()
        
        if options['apply']:
            # Apply MotherPlates
            #-----------------------------------------------
            if hasattr(djTP,'motherplate_ids'):
                print(f" [{self.command_process}] Testplate Apply MotherPlates")
                djTP.clear_cmpbatch_data()
                mp_ids = getattr(djTP,'motherplate_ids')
                for mp in mp_ids:
                    if mp != '-':
                        djMP = MasterPlate.get(mp)
                        if djMP is None:
                            _Valid = False
                            print(f" [{self.command_process}] MasterPlateID {mp} not found")
                        else:                            
                            add_mother_to_testplate(djMP,djTP)      
                djTP.update_n('n_samples')
            
            # Apply Layout
            #-----------------------------------------------
            if _Valid and djTP.n_wells > 0 and djTP.n_reads > 0 and djTP.control_layout:
                print(f" [{self.command_process}] Testplate Apply Layout")
                djTP.n_layout = djTP.apply_layout()            
                
            # Calc Inhibition
            #-----------------------------------------------
            if _Valid and djTP.n_layout > 0:
                print(f" [{self.command_process}] Testplate Calculate Inhibition")
                djTP.calc_inhibition()
        
        if _Valid and options['upload'] and options['overwrite']:
            # Reset AssayData - MIC, CC50, HC50
            #-----------------------------------------------
            if str(djTP.result_type) in DR_CLASSES:
                print(f" [{self.command_process}] Reset AssayData_{str(djTP.result_type)}")
                reset_testplate_assaydata(djTP, upload=options['upload'], overwrite=options['overwrite'],verbose=2)
                
        if _Valid:
            # Process TestPlate -> AssayData
            #-----------------------------------------------
            if str(djTP.result_type) in DR_CLASSES:
                _dr_list = calc_testplate_doseresponse(djTP,upload=options['upload'], overwrite=options['overwrite'],verbose=0)
            elif str(djTP.result_type) in SYN_CLASSES:
                _dr_dict, _iso_lst, _syn_lst = calc_testplate_synergy(djTP,upload=options['upload'], overwrite=options['overwrite'], verbose=0)
            
                
    #------------------------------------------------------------------  
    def handle(self, *args, **options):
    #------------------------------------------------------------------  

        if options['debug']:
            settings.DEBUG = True
        else:
            settings.DEBUG = False
            
        #valLog = Validation_Log("Process_TestPlates")
        
        _dr_list = []
        if options['plateid']:
            djTP = TestPlate.get(options['plateid'])
            self._process_testplate(djTP, *args, **options)
            
            # if djTP is None:
            #     #print(f" [{self.command_process}] PlateID {options['plateid']} not found")
            #     return()
            
            # #logger.info(f" [{self.command_process}] Processing PlateID {djTP} [{djTP.result_type}] ...")
            # if str(djTP.result_type) in DR_CLASSES:
            #     _dr_list = process_testplate_doseresponse(djTP,upload=options['upload'], overwrite=options['overwrite'],verbose=2)
            # elif str(djTP.result_type) in SYN_CLASSES:
            #     _dr_dict, _iso_lst, _syn_lst = process_testplate_synergy(djTP,upload=options['upload'], overwrite=options['overwrite'], verbose=1)
                            
        elif options['runid']:
            djRun = Screen_Run.get(options['runid'])
            if djRun:
                qryTP = TestPlate.objects.filter(run_id=djRun).values_list('plate_id')
                nTP = qryTP.count()
                _iTP = 0
                print(f" [{self.command_process}] Processing RunID {djRun} with {nTP} testplates ...")
                for _tp_id in tqdm(qryTP, total=nTP, desc='TestPlates'):
                    _iTP += 1
                    # djTP = TestPlate.get(_tp_id[0])
                    # #logger.info(f" [{djTP.plate_id}] {djTP.result_type} {djTP.assay_id} {_iTP}/{nTP}")
                    # if str(djTP.result_type) in DR_CLASSES:
                    #     _dr_list = process_testplate_doseresponse(djTP,upload=options['upload'], overwrite=options['overwrite'],verbose=0)
                    # elif str(djTP.result_type) in SYN_CLASSES:
                    #     _dr_dict, _iso_lst, _syn_lst = process_testplate_synergy(djTP,upload=options['upload'], overwrite=options['overwrite'], verbose=1)
            else:
                print(f" [{self.command_process}] RunID {options['runid']} NOT Found")
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

