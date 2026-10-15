import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dsample.models import Project
from applib.django.foreignkey import rename_ForeignKey

class Command(BaseCommand):
    help = "Rename Projects"

    def add_arguments(self, parser):
        parser.add_argument("--oldid")
        parser.add_argument("--newid")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--delete",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Rename_Project")
        
        print(f" [Rename Project] [{options['oldid']}] -> [{options['newid']}] - Upload : {options['upload']} Delete Old: {options['delete']}")
        rename_ForeignKey(Project,options['oldid'],options['newid'],upload=options['upload'],delete=options['delete'])
        
        
