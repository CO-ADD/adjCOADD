import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dplate.models import TestPlate
from applib.bio.doseresponse import process_testplate_doseresponse

class Command(BaseCommand):
    
    command_process = "Process TestPlates"

    def add_arguments(self, parser):
        # parser.add_argument("--projectid")
        # parser.add_argument("--xlsx_file")
        parser.add_argument("--runid",action="store",default=None,help='RunID')
        parser.add_argument("--plateid",action="store",default=None,help='PlateID')
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Process_TestPlates")
        
        _dr_list = []
        if options['plateid']:
            djTP = TestPlate.get(options['plateid'])
            if djTP is None:
                print(f" [{self.command_process}] PlateID {options['plateid']} not found")
                return()
            
            print(f" [{self.command_process}] Processing PlateID {options['plateid']} ...")
            _dr_list = process_testplate_doseresponse(djTP)
            
        elif options['runid']:
            pass
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

