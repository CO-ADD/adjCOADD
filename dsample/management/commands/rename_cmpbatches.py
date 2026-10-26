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
        
    #----------------------------------------------------------------
    def _single_rename(self, oldid,nTotal, nCurrent, upload=False, delete=True):
    #----------------------------------------------------------------
        _p = oldid.split('_')
        if len(_p) == 3:
            newid = ABase_Compound_Batch.str_BatchID("_".join(_p[:2]), _p[2])
            print(f"--- ABase Batch OldID: {oldid} NewID: {newid} --- {nCurrent:8d}/{nTotal:8d}")
            rename_ForeignKey(Compound_Batch,oldid,newid,upload=upload,delete=False)
            #rename_ForeignKey(ABase_Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
            
            # Clean up ABase_Compound_Batch
            if upload:
                _objAbase = ABase_Compound_Batch.get(oldid)
                _objCmp = Compound_Batch.get(oldid)
                    
                if _objAbase:
                    if delete:
                        _objAbase.delete()
                    else:
                        _objAbase.remove()

                if _objCmp:
                    if delete:
                        _objCmp.delete()
                    else:
                        _objCmp.remove()
    
    #----------------------------------------------------------------
    def handle(self, *args, **options):
    #----------------------------------------------------------------
        
        valLog = Validation_Log("Rename_ABase_Batch")
        
        print(f" [Rename ABase Batches] - Upload : {options['upload']} Delete Old: {options['delete']}")
        
        qryABase = ABase_Compound_Batch.objects.all()
        
        nQry = qryABase.count()
        nCmp = 0
        for _obj in qryABase:
            oldid = str(_obj.cmpbatch_id)
            nCmp += 1
            if len(oldid.split('_')) == 3:
                self._single_rename(oldid, nQry,nCmp, upload=options['upload'], delete=options['delete'])
        
        #oldid = 'MCC_006442_011'
        #self._single_rename(oldid, upload=options['upload'], delete=options['delete'])
        
        # _p = oldid.split('_')
        # if len(_p) == 3:
        #     newid = ABase_Compound_Batch.str_BatchID("_".join(_p[:2]), _p[2])
        #     print(f"--- ABase Batch OldID: {oldid} NewID: {newid}")
        #     rename_ForeignKey(Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
        #     #rename_ForeignKey(ABase_Compound_Batch,oldid,newid,upload=options['upload'],delete=options['delete'])
            
        #     # Clean up ABase_Compound_Batch
        #     if options['upload']:
        #         _objAbase = ABase_Compound_Batch.get(oldid)
        #         if _objAbase:
        #             _objAbase.cmpbatch_id = newid
        #             if options['delete']:
        #                 _objAbase.delete()
        #             else:
        #                 _objAbase.remove()

                    
