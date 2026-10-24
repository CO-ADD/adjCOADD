import os
import pandas as pd
import csv
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dsample.models import Project, ABase_Compound_Batch, Compound_Batch
from applib.django.foreignkey import rename_ForeignKey

class Command(BaseCommand):
    help = "Rename ABase Batch ['_' to ':']"

    def add_arguments(self, parser):
        #parser.add_argument("--csvfile")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--delete",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Rename_ABase_Batch")
        
        print(f" [Rename ABase Batches] - Upload : {options['upload']} Delete Old: {options['delete']}")
        
        oldid = 'MCC_000177'
        _p = oldid.split('_')
        if len(_p) == 3:
            newid = ABase_Compound_Batch.str_BatchID("_".join(_p[:2]), _p[2])
            print(f"ABase Batch OldID: {oldid} NewID: {newid}")
            rename_ForeignKey(ABase_Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
            rename_ForeignKey(Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])

        # _abase_batch = ABase_Compound_Batch.objects.all()     
        # for _b in _abase_batch:
        #     oldid = str(_b.cmpbatch_id)
        #     _p = oldid.split('_')
        #     if len(_p) == 3:
        #         newid = ABase_Compound_Batch.str_BatchID("_".join(_p[:2]), _p[2])
        #         print(f"ABase Batch OldID: {oldid} NewID: {newid}")
        #         rename_ForeignKey(ABase_Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
        #         rename_ForeignKey(Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
                    
