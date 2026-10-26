
from django.apps import apps

from django.db.models import ForeignKey, Model
from django.db.models.base import ModelBase
from applib.data.set_fielddata import rename_fkeyarray
from tqdm import tqdm

#-----------------------------------------------------------------------------------
def get_Models_byForeignKey(fkModel):
#-----------------------------------------------------------------------------------

    FOREIGNKEY_ARRAY_FIELDS = ['cmpbatch_lst']

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
                if f.name == 'cmpbatch_id' and hasattr(m,'cmpbatch_lst'):
                    models_with_fk.append({'Model':m,'Field':f.name, 'ArrayField':'cmpbatch_lst'})
                else:
                    models_with_fk.append({'Model':m,'Field':f.name, 'ArrayField':None})    

    return(models_with_fk)

#-----------------------------------------------------------------------------------
def update_ForeignKeyList(fklModel, oldPK, newPK, use_temp_pk=False, upload=False):
#-----------------------------------------------------------------------------------
    print(f" [update_ForeignKeyList] {fklModel.__name__} {oldPK} -> {newPK} [Upload:{upload}]")
    
    
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
                
            for _fkModel in fkModel_Lst:
                
                # For each fkModel get objects with foreignkey = oldPK
                filter_params = {_fkModel['Field']: oldPK}
                qryFK = _fkModel['Model'].objects.filter(**filter_params)
                _nqry = qryFK.count()
                
                _ostr = f" [rename_ForeignKey] {_fkModel['Model'].__name__}.{_fkModel['Field']} -> {_nqry}"
                if _nqry > 0:
                    for fkObj in tqdm(qryFK, total=_nqry, desc=_ostr):
                        setattr(fkObj,_fkModel['Field'],djNew)
                        
                        if upload:
                            fkObj.save()
                else:
                    print(_ostr)
                    
                # Update ForeignKeyArray Fields
                if _fkModel['ArrayField'] is not None:
                    filter_params = {_fkModel['ArrayField']+'__contains':[oldPK]}
                    qryFK = _fkModel['Model'].objects.filter(**filter_params)
                    _nqry = qryFK.count()
                    
                    _ostr = f" [rename_ForeignKey] {_fkModel['Model'].__name__}.{_fkModel['ArrayField']} -> {_nqry}"
                    if _nqry > 0:
                        for fkObj in tqdm(qryFK, total=_nqry, desc=_ostr):
                            rename_fkeyarray(fkObj, oldPK, newPK, _fkModel['ArrayField'])
                            
                            if upload:
                                fkObj.save()
                    else:
                        print(_ostr)
                        
                        
                
                        
            # Remove/Delete OldPK
            if upload:
                if delete:
                    djOld.delete()
                else:
                    djOld.remove()
                    
        else:
            print(f' [rename_ForeignKey] Warning: NewPK {newPK} Exists ')
    else:
        print(f' [rename_ForeignKey] Error: OldPK {oldPK} Dose NOT Exists')

