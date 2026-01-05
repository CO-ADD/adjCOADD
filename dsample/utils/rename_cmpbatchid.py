from django.apps import apps

from django.db.models import ForeignKey, Model
from django.db.models.base import ModelBase
from applib.django.foreignkey import get_Models_byForeignKey

from dsample.models import Compound_Batch
from dscreen.models import AssayData_CC50, AssayData_HC50, AssayData_MIC
from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Inhib, Summary_CmpBatch_Doseresp
from dplate.models import TestWell, MasterWell

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------
def rename_CmpBatchID(oldCmpBatchID, newCmpBatchID=None, upload=False, remove=True,):
#-----------------------------------------------------------------------------------

    CMPBATCHLST_MODELS = [AssayData_CC50, AssayData_HC50, AssayData_MIC, 
                          TestWell, MasterWell, 
                          Summary_CmpBatch, Summary_CmpBatch_Inhib, Summary_CmpBatch_Doseresp]

    print(f" [rename_CmpBatchID]  {oldCmpBatchID} -> {newCmpBatchID} [Upload:{upload} Remove:{remove}]")
    djOld = Compound_Batch.get(oldCmpBatchID)
    djNew = Compound_Batch.get(newCmpBatchID)

    if djOld:
        if djNew is None:
            # Get list of Models with fkModel as ForeignKey
            fkModel_Lst = get_Models_byForeignKey(Compound_Batch)
            
            # Create NewPK
            djNew = Compound_Batch.get(oldCmpBatchID)
            djNew.pk = newCmpBatchID
            if upload:
                djNew.save()
                
            for fkModel in fkModel_Lst:
                # For each fkModel get objects with foreignkey = djOld
                filter_params = {fkModel['Field']: djOld}
                qryFK = fkModel['Model'].objects.filter(**filter_params)
                print(f" [rename_CmpBatchID] { fkModel['Model'].__name__} -> {qryFK.count()} ") 
                for fkObj in qryFK:
                    setattr(fkObj,fkModel['Field'],djNew)
                    if upload:
                        fkObj.save()

            # Remove/Delete OldPK
            if upload:
                if remove:
                    djOld.remove()
                else:
                    djOld.delete()

            # Update CmpBatch_Lst
            for lstModel in CMPBATCHLST_MODELS:
                qryCmpLst = lstModel.objects.filter(cmpbatch_lst__contains=[oldCmpBatchID])
                print(f" [rename_CmpBatchID] Lst { lstModel['Model'].__name__} -> {qryCmpLst.count()} ") 
                for djInst in qryCmpLst:
                    newCmpBatchLst = []
                    newLst = False
                    for cmpbatch in djInst.cmpbatch_lst:
                        if cmpbatch == oldCmpBatchID:
                            newCmpBatchLst.append(newCmpBatchID)
                            newLst = True
                        else:
                            newCmpBatchLst.append(oldCmpBatchID)
                    djInst.cmpbatch_lst = newCmpBatchLst

                    if newLst and upload:
                        djInst.save()

        else:
            print(f' [rename_OrgBatchID] Error: New CmpBatchID {newCmpBatchID} Exists')
    else:
        print(f' [rename_OrgBatchID] Error: Old CmpBatchID {oldCmpBatchID} Dose NOT Exists')


#-----------------------------------------------------------------------------------
def rename_cmpbatch_lst():
#-----------------------------------------------------------------------------------
    pass