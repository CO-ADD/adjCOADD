from django.apps import apps

from django.db.models import ForeignKey, Model
from django.db.models.base import ModelBase
from applib.django.foreignkey import get_Models_byForeignKey

from dorganism.models import Organism_Batch
from dorganism.utils.utils import get_subdir

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------
def rename_OrgBatchID(oldOrgBatchID, newOrgBatchID=None, upload=False, remove=True,):
#-----------------------------------------------------------------------------------
    print(f" [rename_OrgBatchID]  {oldOrgBatchID} -> {newOrgBatchID} [Upload:{upload} Remove:{remove}]")
    djOld = Organism_Batch.get(oldOrgBatchID)
    djNew = Organism_Batch.get(newOrgBatchID)
    
    if djOld:
        if djNew is None:
            # Get list of Models with fkModel as ForeignKey
            fkModel_Lst = get_Models_byForeignKey(Organism_Batch)
            
            # Create NewPK
            djNew = Organism_Batch.get(oldOrgBatchID)
            djNew.pk = newOrgBatchID
            if upload:
                djNew.save()
                
            for fkModel in fkModel_Lst:
                # For each fkModel get objects with foreignkey = djOld
                filter_params = {fkModel['Field']: djOld}
                qryFK = fkModel['Model'].objects.filter(**filter_params)
                print(f" [rename_OrgBatchID] { fkModel['Model'].__name__} -> {qryFK.count()} ") 
                for fkObj in qryFK:
                    setattr(fkObj,fkModel['Field'],djNew)


                    if fkModel.__name__ == 'OrgBatch_Image':
                        # Rename orgbatch_ID Images
                        _oldI = fkObj.image_name
                        _old = fkObj.image_file
                        fkObj.image_name = _oldI.replace(oldOrgBatchID,newOrgBatchID)
                        fkObj.image_file = f"images/orgbatch/{get_subdir(fkObj.image_name)}/{fkObj.image_name}"
                        logger.warning(f"[  XXX Rename {fkModel.__name__} --> Rename Image Files {_old} -> {fkObj.image_file}")
                        
                    elif fkModel.__name__ == 'Genome_Sequence':
                        # Rename orgbatch_ID Sequences
                        _old = fkObj.seq_name
                        fkObj.seq_name = _old.replace(oldOrgBatchID,newOrgBatchID)
                        fkObj.source_code = fkObj.source_code.replace(oldOrgBatchID,newOrgBatchID)
                        fkObj.source_link = fkObj.source_link.replace(oldOrgBatchID,newOrgBatchID)
                        logger.warning(f"[  XXX Rename {fkModel.__name__} --> Rename Sequence Folder/Files {_old} -> {fkObj.seq_name}")

                    if upload:
                        fkObj.save()

            # Remove/Delete OldPK
            if upload:
                if remove:
                    djOld.remove()
                else:
                    djOld.delete()
                    
        else:
            print(f' [rename_OrgBatchID] Error: New OrgBatchID {newOrgBatchID} Exists')
    else:
        print(f' [rename_OrgBatchID] Error: Old OrgBatchID {oldOrgBatchID} Dose NOT Exists')

