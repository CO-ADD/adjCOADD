import os, sys
import pandas as pd
import logging

from decimal import Decimal

from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dplate.models import TestPlate
from dscreen.models import Screen_Run
from adjcoadd.constants import DR_CLASSES, SYN_CLASSES
from applib.bio.doseresponse import process_testplate_doseresponse

# from applib.bio.bio_data import DR_Range
# from applib.bio.synergy import Synergy_TestPlate, Synergy
# from applib.bio.doseresponse import DoseResponse
# from dplate.models import TestPlate, TestWell
# from dscreen.models import AssayData_MIC, AssayData_SynergyIsobol 
from adjcoadd.constants import SYN_CLASSES, COMPOUND_SEP
from applib.bio.synergy import process_testplate_synergy

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
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        
    #------------------------------------------------------------------  
    def handle(self, *args, **options):
    #------------------------------------------------------------------  

        # Set Loggers and Validation_Log        
        valLog = Validation_Log("Process_TestPlates")
        # if hasattr(logging, 'getHandlerByName'):
        #     console_handler = logging.getHandlerByName('console')
        # else:
        #     console_handler = logging.StreamHandler(sys.stdout)
        # if console_handler not in logger.handlers:
        #     logger.addHandler(console_handler)
        # logger.setLevel(logging.DEBUG)
        #------------------------------------------------
        
        _dr_list = []
        if options['plateid']:
            djTP = TestPlate.get(options['plateid'])
            if djTP is None:
                #print(f" [{self.command_process}] PlateID {options['plateid']} not found")
                return()
            
            #logger.info(f" [{self.command_process}] Processing PlateID {djTP} [{djTP.result_type}] ...")
            if str(djTP.result_type) in DR_CLASSES:
                _dr_list = process_testplate_doseresponse(djTP,verbose=2)
            elif str(djTP.result_type) in SYN_CLASSES:
                _dr_dict, _iso_lst, _syn_lst = process_testplate_synergy(djTP,upload=options['upload'], overwrite=options['overwrite'], verbose=2)
                            
        elif options['runid']:
            djRun = Screen_Run.get(options['runid'])
            if djRun:
                qryTP = TestPlate.objects.filter(run_id=djRun).values_list('plate_id')
                nTP = qryTP.count()
                _iTP = 0
                logger.info(f" [{self.command_process}] Processing RunID {djRun} with {nTP} testplates ...")
                for _tp_id in qryTP:
                    _iTP += 1
                    djTP = TestPlate.get(_tp_id[0])
                    logger.info(f" [{djTP.plate_id}] {djTP.result_type} {djTP.assay_id} {_iTP}/{nTP}")
                    if str(djTP.result_type) in DR_CLASSES:
                        _dr_list = process_testplate_doseresponse(djTP,verbose=0)
                    elif str(djTP.result_type) in SYN_CLASSES:
                        _dr_dict, _iso_lst, _syn_lst = process_testplate_synergy(djTP,upload=options['upload'], overwrite=options['overwrite'], verbose=2)
            
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

