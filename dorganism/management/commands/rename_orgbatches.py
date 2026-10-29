import os
import pandas as pd
import csv
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dorganism.models import Organism_Batch
from applib.django.foreignkey import rename_ForeignKey

class Command(BaseCommand):
    help = "Rename Organism Batch ['_' to ':']"

    def add_arguments(self, parser):
        #parser.add_argument("--csvfile")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--delete",action="store_true",default=False)
        
    #----------------------------------------------------------------
    def _single_rename(self, oldid, nTotal, nCurrent, upload=False, delete=True):
    #----------------------------------------------------------------
        # 001 - Changing from MCC_NNNNNN_BBB to MCC_NNNNNN:BBB
        _p = oldid.split('_')
        if len(_p) > 2:
            # 001
            newid = Organism_Batch.str_OrgBatchID("_".join(_p[:2]), _p[2])
        
            print(f"--- Org Batch OldID: {oldid} NewID: {newid} --- {nCurrent:8d}/{nTotal:8d}")
            rename_ForeignKey(Organism_Batch,oldid,newid,upload=upload,delete=False)
    
    #----------------------------------------------------------------
    def handle(self, *args, **options):
    #----------------------------------------------------------------
        
        valLog = Validation_Log("Rename_Organism_Batch")
        
        print(f" [Rename OrgBatch Batches] - Upload : {options['upload']} Delete Old: {options['delete']}")
        
        qryOrg = Organism_Batch.objects.all()
        
        nQry = qryOrg.count()
        nBatch = 0
        for _obj in qryOrg:
            oldid = str(_obj.orgbatch_id)
            nBatch += 1
            if len(oldid.split('_')) > 1:
                self._single_rename(oldid, nQry,nBatch, upload=options['upload'], delete=options['delete'])
        
                    
