import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dplate.models import TestPlate


class Command(BaseCommand):
    
    command_process = "Process TestPlates"

    def add_arguments(self, parser):
        # parser.add_argument("--projectid")
        # parser.add_argument("--xlsx_file")
        parser.add_argument("--runid",action="store",default=None,help='RunID')
        parser.add_argument("--plateid",action="store",default=None,help='PlateID')
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Process_TestPlates")
        
        if options['plateid']:
            pass
            
        elif options['runid']:
            pass
        else:
            print(f" [{self.command_process}] --plateid <PlateID> or --runid <RunID> ")
            print(self.help)
            return()
        
        

