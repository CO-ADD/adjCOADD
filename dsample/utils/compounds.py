from dplate.models import MasterPlate, TestPlate, TestWell
from dsample.models import Project, COADD_Compound, Compound_Batch, ABase_Compound_Batch

from applib.data.set_fielddata import set_model_arrayfields, set_model_fields, set_model_dicts

#-----------------------------------------------------------------------------------
def Upload_COADD_Compound(CompoundDict, upload=False, overwrite=False, appuser=None, valLog=None, verbose=0):
#-----------------------------------------------------------------------------------

    cpyFields = ['compound_code','compound_name','compound_description',
                    'cpoz_sn','cpoz_id','coadd_id','spark_id',
                'reg_smiles','reg_structure','reg_mw','reg_mf',
                'reg_amount','reg_volume','reg_conc','reg_solvent',
                'stock_amount','prep_date',
                'pub_status','pub_date',
                'ora_compound_id','ora_project_id','ora_compound_type'
                ]

    dictFields = [
                'reg_amount_unit','reg_volume_unit','reg_conc_unit','stock_amount_unit',
                ]
    
    # - COADD_Compound ----------------------        
    if 'compound_id' in CompoundDict:        
        djCmpd = COADD_Compound.get(CompoundDict['compound_id'])
    else:
        djCmpd = None

    if djCmpd is None:
        djCmpd = COADD_Compound()
        new_compound = True
        sample_status = 'New'
        if 'compound_id' in CompoundDict:
            djCmpd.compound_id = CompoundDict['compound_id']
    else:
        new_compound = False
        sample_status = 'ID exists'
        if valLog:
            valLog.add_warning()

        #outNumbers['New Compounds'] += 1
    #else:
        #row['Issue'] = f"Exists"

    djPrj= Project.get(CompoundDict['project_id'])
    djCmpd.project_id = djPrj

    set_model_fields(djCmpd,CompoundDict,cpyFields)
    set_model_dicts(djCmpd,CompoundDict,dictFields)

    if djCmpd.reg_mw < 2:
        djCmpd.reg_mw = 0
    if djCmpd.reg_mf == 'CxHxNxOx':
        djCmpd.reg_mf = ''    

    # - Validation ----------------------------------------------------------
    validStatus = True
    djCmpd.set_defaults_model()
    validDict = djCmpd.validate_fields()
    if validDict:
        validStatus = False
        for k in validDict:
            print('Warning',k,validDict[k],'-')
        #outDict.append(row)

    # Check if compound_code already exists in djPrj
    qryCode = COADD_Compound.objects.filter(compound_code=djCmpd.compound_code, project_id=djPrj)
    if qryCode.count() > 0:
        sample_status = 'Code exists'
        validStatus = False
        valLog.add_error('Duplicate Code',djCmpd.compound_code,'Code exists already','Correct Compound Code')

    print(f"[{djCmpd.compound_code}] - {djCmpd.reg_smiles} {djCmpd.reg_mw} - {djPrj} - {sample_status}")

    if validStatus:
        if upload:
            if new_compound or overwrite:
                #print(f"[{djCmpd.compound_code}] Uploading")
                #outNumbers['Upload Batches'] += 1
                djCmpd.save()

 