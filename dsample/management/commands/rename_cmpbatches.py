import os
import pandas as pd
import csv
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dsample.models import Project
from applib.django.foreignkey import rename_ForeignKey

class Command(BaseCommand):
    help = "Rename Projects"

    def add_arguments(self, parser):
        parser.add_argument("--csvfile")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--delete",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Rename_Project")
        
        print(f" [Rename CmpBatches] [{options['csvfile']}] - Upload : {options['upload']} Delete Old: {options['delete']}")
        
        if options['csvfile']:
            if os.path.isfile(options['csvfile']):                
                with open(options['csvfile'], mode='r') as file:
                    reader = csv.DictReader(file)
                    for row in reader:
                        rename_ForeignKey(Project,row['oldid'],row['newid'],upload=options['upload'],delete=options['delete'])