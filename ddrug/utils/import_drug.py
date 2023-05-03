import os
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST, MIC_COADD, MIC_Pub
from ddrug.utils.molecules import *

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import *

# ----------------------------------------------------------------------------------------------------
def imp_Drug_fromDict(iDict,valLog):
    """
    Create Drug instance from {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    # Find Instance if exist
    djDrug = Drug.get(iDict['drug_name'],None)
    if djDrug is None:
        djDrug = Drug()
        #djDrug.drug_id = iDict['drug_id']
        djDrug.drug_name = iDict['drug_name']
        valLog.add_log('Info','New Drug',f"{iDict['drug_name']}",'-')
    djDrug.drug_type = Dictionary.get(djDrug.Choice_Dictionary["drug_type"],iDict['drug_type'])

    djDrug.n_compounds = iDict['ncmpd']
    djDrug.drug_othernames = split_StrList(iDict['drug_othernames'])
    djDrug.drug_codes = split_StrList(iDict['drug_code'])
    djDrug.drug_note = iDict['drug_note']
    djDrug.drug_panel = split_StrList(iDict['panel'])

    djDrug.antimicro = iDict['antimicro']
    djDrug.antimicro_class = iDict['antimicro_class']
    djDrug.drug_class = iDict['drug_class']
    djDrug.drug_subclass = iDict['drug_subclass']
    djDrug.drug_target = iDict['drug_target']
    #djDrug.drug_subtarget = iDict['drug_subtarget']
    #djDrug.moa = iDict['moa']

    djDrug.vendor = iDict['vendor']
    djDrug.vendor_catno = iDict['vendor_catno']

    djDrug.chembl = iDict['chembl']
    djDrug.drugbank = iDict['drugbank']
    djDrug.cas = iDict['cas']
    djDrug.pubchem = iDict['pubchem']
    djDrug.chemspider = iDict['chemspider']
    djDrug.unii = iDict['unii']
    djDrug.kegg = iDict['kegg']
    djDrug.comptox = iDict['comptox']
    djDrug.echa = iDict['echa']
    djDrug.chebi = iDict['chebi']
    djDrug.uq_imb = iDict['imb']

    djDrug.smiles = iDict['smiles']
    djDrug.smol = iDict['smiles']
    djDrug.smol = smiles2mol(iDict['smiles'],verbose=1)

    djDrug.clean_Fields()
    validStatus = True
    validDict = djDrug.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djDrug.VALID_STATUS = validStatus

    return(djDrug)


# ----------------------------------------------------------------------------------------------------
def imp_VitekCard_fromDict(iDict,valLog):
    """
    Create VITEK_Card instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    
    djVitekCard = VITEK_Card.get(iDict['card_barcode'])
    if djVitekCard is None:
        djVitekCard = VITEK_Card()
        djVitekCard.card_barcode = iDict['card_barcode']
        valLog.add_log('Info','New VITEK card',f"{iDict['card_barcode']}-{iDict['card_code']}",'-')
    else:
        valLog.add_log('Info','Update VITEK card',f"{iDict['card_barcode']} -{iDict['card_code']}",'-')

    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error','Organism Batch does not Exists',iDict['orgbatch_id'],'Use existing OrganismBatch ID')
        validStatus = False
    djVitekCard.orgbatch_id = OrgBatch

    djVitekCard.card_type = Dictionary.get(djVitekCard.Choice_Dictionary["card_type"],iDict['card_type'])
    if djVitekCard.card_type is None:
        valLog.add_log('Error','Vitek Card Type not Correct',iDict['card_type'],'-')
        validStatus = False

    djVitekCard.card_code = iDict['card_code']
    djVitekCard.instrument = iDict['instrument']
    djVitekCard.expiry_date = iDict['expiry_date']
    djVitekCard.proc_date = iDict['processing_date']
    djVitekCard.analysis_time = iDict['analysis_time']

    djVitekCard.clean_Fields()
    validDict = djVitekCard.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djVitekCard.VALID_STATUS = validStatus
    return(djVitekCard)

# ----------------------------------------------------------------------------------------------------
def imp_VitekID_fromDict(iDict,valLog):
    """
    Create VITEK_ID instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    Barcode = VITEK_Card.get(iDict['card_barcode']) 
    if Barcode is None:
        validStatus = False
        valLog.add_log('Error','VITEK card does not Exists',f"{iDict['card_code']} ({iDict['card_barcode']})",'-')

    djVitekID = VITEK_ID.get(Barcode)
    if djVitekID is None:
        djVitekID = VITEK_ID()
        djVitekID.card_barcode = Barcode
        valLog.add_log('Info','New VITEK ID',f"{iDict['card_code']} ({iDict['card_barcode']})",'-')
    else:
        valLog.add_log('Info','Update VITEK ID',f"{iDict['card_code']} ({Barcode})",'-')

    djVitekID.process = iDict['vitek_process']
    djVitekID.id_organism = iDict['id_organism']
    djVitekID.id_probability = iDict['id_probability']
    djVitekID.id_confidence = iDict['id_confidence']
    #retInstance.id_source = iDict['card_barcode']
    djVitekID.filename = iDict['filename']
    djVitekID.page_no = iDict['pageno']  

    djVitekID.clean_Fields()
    validDict = djVitekID.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')


    djVitekID.VALID_STATUS = validStatus
    return(djVitekID)


# ----------------------------------------------------------------------------------------------------
def imp_VitekAST_fromDict(iDict,valLog):
    """
    Create VITEK_ID instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    Barcode = VITEK_Card.get(iDict['card_barcode']) 
    if Barcode is None:
        validStatus = False
        valLog.add_log('Error','VITEK card does not Exists',f"{iDict['card_code']} ({iDict['card_barcode']})",'-')

    DrugID = Drug.get(iDict['drug_name'])
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error','Drug does not Exists',f"{iDict['drug_name']} ({iDict['card_barcode']})",'-')

    if validStatus:
        djVitekAST = VITEK_AST.get(Barcode,DrugID,iDict['bp_source'],iDict['selected_organism'])
    else:
        djVitekAST = None
            
    if djVitekAST is None:
        djVitekAST = VITEK_AST()
        djVitekAST.card_barcode = Barcode
        djVitekAST.drug_id = DrugID
        djVitekAST.bp_source = iDict['bp_source']
        valLog.add_log('Info','New VITEK AST',f"{Barcode} {DrugID} {iDict['bp_source']}",'-')
    
    djVitekAST.mic = iDict['mic']
    djVitekAST.process = iDict['vitek_process']
    djVitekAST.bp_profile = iDict['bp_profile']
    djVitekAST.bp_comment = iDict['bp_comment']
    djVitekAST.selection = iDict['organism_origin']
    djVitekAST.organism = iDict['selected_organism']
    djVitekAST.filename = iDict['filename']
    djVitekAST.page_no = iDict['pageno']  

    djVitekAST.clean_Fields()
    validDict = djVitekAST.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djVitekAST.VALID_STATUS = validStatus
    return(djVitekAST)
 
# ----------------------------------------------------------------------------------------------------
def imp_MICCOADD_fromDict(iDict,valLog):
    """
    Create VITEK_ID instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    DrugID = Drug.get(iDict['drug_name'])
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error','Drug does not Exists',f"{iDict['drug_name']} ",'-')

    OrgBatchID = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatchID is None:
        validStatus = False
        valLog.add_log('Error','OrgBatchID does not Exists',f"{iDict['orgbatch_id']} ",'-')

    if validStatus:
        djMIC = MIC_COADD.get(OrgBatchID,DrugID,iDict['testplate_id'],iDict['testwell_id'])
    else:
        djMIC = None
            
    if djMIC is None:
        djMIC = MIC_COADD()
        djMIC.orgbatch_id = OrgBatchID
        djMIC.drug_id = DrugID
        djMIC.run_id = iDict['run_id']
        djMIC.testplate_id = iDict['testplate_id']
        djMIC.testwell_id = iDict['testwell_id']
        valLog.add_log('Info','New MIC ',f"{OrgBatchID} {DrugID} {iDict['testplate_id']}:{iDict['testwell_id']}",'-')
    
    djMIC.mic = iDict['mic']
    djMIC.mic_unit = iDict['mic_unit']
    djMIC.mic_type = Dictionary.get(MIC_COADD.Choice_Dictionary["mic_type"],'BMD',None,verbose=1)

    djMIC.plate_size = Dictionary.get(MIC_COADD.Choice_Dictionary["plate_size"],iDict['plate_size'],None,verbose=1)
    djMIC.plate_material = Dictionary.get(MIC_COADD.Choice_Dictionary["plate_material"],iDict['plate_material'],None,verbose=1)

    #djMIC.bp_profile = iDict['bp_profile']
    #djMIC.bp_source = iDict['bp_source']
    #djMIC.media = Dictionary.get(cls.Choice_Dictionary["media"],iDict['media'],None,verbose=1)

    djMIC.clean_Fields()
    validDict = djMIC.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djMIC.VALID_STATUS = validStatus
    
    return(djMIC)

# ----------------------------------------------------------------------------------------------------
def imp_MICPub_fromDict(iDict,valLog):
    """
    Create MIC_Pub instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    DrugID = Drug.get(iDict['drug_name'].strip())
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error','Drug does not Exists',f"{iDict['drug_name']} ",'-')

    OrganismID = Organism.get(iDict['organism_id']) 
    if OrganismID is None:
        validStatus = False
        valLog.add_log('Error','OrganismID does not Exists',f"{iDict['organism_id']} ",'-')

    if validStatus:
        djMIC = MIC_Pub.get(OrganismID,DrugID,iDict['source'])
    else:
        djMIC = None
            
    if djMIC is None:
        djMIC = MIC_Pub()
        djMIC.organism_id = OrganismID
        djMIC.drug_id = DrugID
        djMIC.source = iDict['source']
        valLog.add_log('Info','New MIC ',f"{iDict['organism_id']} {iDict['drug_name']} {iDict['source']}",'-')
    
    if 'source_type' in iDict:
        djMIC.mic_type = Dictionary.get(MIC_Pub.Choice_Dictionary["mic_type"],iDict['source_type'],None,verbose=1)
    if 'mic' in iDict:
        djMIC.mic = iDict['mic']
        djMIC.mic_unit = iDict['mic_unit']
    if 'bp_profile' in iDict:
        djMIC.bp_profile = iDict['bp_profile']
        if 'bp_source' in iDict:
            djMIC.bp_source = iDict['bp_source']
    if 'zone_diameter' in iDict:
        djMIC.zone_diameter = iDict['zone_diameter']


    djMIC.clean_Fields()
    validDict = djMIC.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djMIC.VALID_STATUS = validStatus
    
    return(djMIC)
