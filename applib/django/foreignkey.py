
from django.apps import apps

from django.db.models import ForeignKey, Model
from django.db.models.base import ModelBase

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
def rename_ForeignKey(fkModel, oldPK, newPK, use_temp_pk=False):
#-----------------------------------------------------------------------------------
    print(f" {oldPK} -> {newPK}")
    djOld = fkModel.get(oldPK)
    if djOld:
        fkmodel_lst = get_Models_byForeignKey(fkModel)
        for fk in fkmodel_lst:
            filter_params = {fk['Field']: oldPK}
            qryFK = fk['Model'].objects.filter(**filter_params)
            print(f" {fk['Model']} {qryFK.count()}")
