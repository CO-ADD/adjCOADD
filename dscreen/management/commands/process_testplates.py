import os, sys
import pandas as pd
import logging

from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dplate.models import TestPlate
from adjcoadd.constants import DR_CLASSES, SYN_CLASSES
from applib.bio.doseresponse import process_testplate_doseresponse
from applib.bio.synergy import process_testplate_synergy

logger = logging.getLogger(__name__)

class Command(BaseCommand):
    
    command_process = "Process TestPlates"

    def add_arguments(self, parser):
        # parser.add_argument("--projectid")
        # parser.add_argument("--xlsx_file")
        parser.add_argument("--runid",action="store",default=None,help='RunID')
        parser.add_argument("--plateid",action="store",default=None,help='PlateID')
        
    def handle(self, *args, **options):

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
                print(f" [{self.command_process}] PlateID {options['plateid']} not found")
                return()
            
            logger.info(f" [{self.command_process}] Processing PlateID {djTP} [{djTP.result_type}] ...")
            if str(djTP.result_type) in DR_CLASSES:
                _dr_list = process_testplate_doseresponse(djTP,verbose=2)
            elif str(djTP.result_type) in SYN_CLASSES:
                _dr_list = process_testplate_synergy(djTP,verbose=2)
                
                print(_dr_list)
            
        elif options['runid']:
            pass
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

