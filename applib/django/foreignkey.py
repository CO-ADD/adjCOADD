
from django.apps import apps

from django.db.models import ForeignKey, Model
from django.db.models.base import ModelBase
from tqdm import tqdm

#-----------------------------------------------------------------------------------
def get_Models_byForeignKey(fkModel):
#-----------------------------------------------------------------------------------
    all_models = apps.get_models()
    #print(all_models)
    #print(f" {fkModel} {type(fkModel)} {Model}")

    objModel = None
    if isinstance(fkModel,(Model, ModelBase)):
        objModel = fkModel

    if isinstance(fkModel,str):
        for m in all_models:
            if m.__name__ == fkModel:
                objModel = m

    # Find models with <fkModel> as a foreign key
    models_with_fk = []
    for m in all_models:
        #print(m.__name__)
        for f in m._meta.fields:
            #print(type(f), f.related_model)
            if isinstance(f, ForeignKey) and f.related_model == objModel:
                models_with_fk.append({'Model':m,'Field':f.name})    

    return(models_with_fk)


#-----------------------------------------------------------------------------------
def rename_ForeignKey(fkModel, oldPK, newPK, use_temp_pk=False, upload=False, delete=True):
#-----------------------------------------------------------------------------------
    print(f" [rename_ForeignKey] {fkModel.__name__} {oldPK} -> {newPK} [Upload:{upload} Remove:{delete}]")
    djOld = fkModel.get(oldPK)
    djNew = fkModel.get(newPK)
    
    if djOld:
        if djNew is None:
            # Get list of Models with fkModel as ForeignKey
            fkModel_Lst = get_Models_byForeignKey(fkModel)
            
            # Create NewPK
            djNew = fkModel.get(oldPK)
            djNew.pk = newPK
            if upload:
                djNew.save()
                
            for fkModel in fkModel_Lst:
                # For each fkModel get objects with foreignkey = oldPK
                filter_params = {fkModel['Field']: oldPK}
                qryFK = fkModel['Model'].objects.filter(**filter_params)
                _nqry = qryFK.count()
                
                print(f" [rename_ForeignKey] { fkModel['Model'].__name__} -> {_nqry} ") 
                for fkObj in tqdm(qryFK, total=_nqry, desc=fkModel['Model'].__name__):
                    setattr(fkObj,fkModel['Field'],djNew)
                    if upload:
                        fkObj.save()
                        
            # Remove/Delete OldPK
            if upload:
                if delete:
                    djOld.delete()
                else:
                    djOld.remove()
                    
        else:
            print(f' [rename_ForeignKey] Error: NewPK {newPK} Exists')
    else:
        print(f' [rename_ForeignKey] Error: OldPK {oldPK} Dose NOT Exists')

