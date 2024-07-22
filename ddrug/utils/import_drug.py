import os
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST, MIC_COADD, MIC_Pub, Breakpoint
from dscreen.models import Screen_Run
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
        valLog.add_log('Info',"",f"{iDict['drug_name']}",'New Drug','-') 
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
#    djDrug.smol = iDict['smiles']
#    djDrug.smol = smiles2mol(iDict['smiles'],verbose=1)

    djDrug.clean_Fields()
    validStatus = True
    validDict = djDrug.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning','',k,validDict[k],'-')

    djDrug.VALID_STATUS = validStatus

    return(djDrug)

# ----------------------------------------------------------------------------------------------------
def imp_Breakpoint_fromDict(iDict,valLog,upload=False):
    """
    Create Breakpoint instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True

    DrugID = Drug.get(iDict['drug_name'])
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error','oraOrgDB',f"{iDict['drug_name']} ",'BP Drug does not Exists','-')


    if 'org_name' in iDict:
        OrgName = iDict['org_name']
        OrgRank = Dictionary.get(Breakpoint.Choice_Dictionary["org_rank"],iDict['org_rank'])
        if OrgRank is None:
            valLog.add_log('Error','oraOrgDB',iDict['org_rank'],'Tax Rank not correct','-')
            validStatus = False
    else:
        OrgName = None
        OrgRank = None

    if 'notorg_name' in iDict:
        NotOrgName = iDict['notorg_name']
        NotOrgRank = Dictionary.get(Breakpoint.Choice_Dictionary["notorg_rank"],iDict['notorg_rank'])
        if NotOrgRank is None:
            valLog.add_log('Error','oraOrgDB',iDict['notorg_rank'],'(Not) Tax Rank not correct','-')
            validStatus = False
    else:
        NotOrgName = None
        NotOrgRank = None

    djBP = Breakpoint.get(DrugID, OrgName, OrgRank, NotOrgName, NotOrgRank,
                        iDict['medical_application'], iDict['bp_type'], iDict['bp_source'])
    if djBP is None:
        djBP = Breakpoint()
        djBP.drug_id = DrugID
        djBP.org_name = OrgName
        djBP.org_rank = OrgRank
        djBP.notorg_name = NotOrgName
        djBP.notorg_rank = NotOrgRank

        valLog.add_log('Info',"",f"{iDict['drug_name']} {OrgRank} {OrgName} {NotOrgRank} {NotOrgName}",'New BP','-')

    djBP.bp_type = Dictionary.get(djBP.Choice_Dictionary["bp_type"],iDict['bp_type'])
    if djBP.bp_type is None:
        valLog.add_log('Error','oraOrgDB',iDict['bp_type'],'BP Type not correct','-')
        validStatus = False

    djBP.med_application = iDict['medical_application']
    djBP.bp_res_gt = iDict['bp_resistant_gt']
    djBP.bp_sens_le = iDict['bp_sensitive_le']
    djBP.bp_unit = iDict['bp_unit']
    djBP.bp_comb = iDict['combination_type']
    djBP.bp_source = iDict['bp_source']
    djBP.bp_source_version = iDict['bp_source_version']

    djBP.clean_Fields()
    validStatus = True
    validDict = djBP.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning','',k,validDict[k],'-')
            #print(f"Warning : {k} {validDict[k]}")
    djBP.VALID_STATUS = validStatus

    return(djBP)


# ----------------------------------------------------------------------------------------------------
def imp_VitekCard_fromDict(iDict,valLog,upload=False):
    """
    Create VITEK_Card instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    
    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error',iDict['filename'],iDict['orgbatch_id'],'Organism Batch does not exists','Use existing OrganismBatch ID')
        validStatus = False

    infoCard = f"Card: {iDict['card_code']} ({iDict['card_barcode']})\n OrgBatch: {iDict['orgbatch_id']}"
    djVitekCard = VITEK_Card.get(iDict['card_barcode'])
    if djVitekCard is None:
        djVitekCard = VITEK_Card()
        djVitekCard.card_barcode = iDict['card_barcode']
        if upload:
            #print(f"New [{upload}] -> Info")
            valLog.add_log('Info',iDict['filename'],infoCard, f"New {iDict['card_type']} VITEK card",'-')
        else:
            #print(f"New [{upload}] -> Warning")
            valLog.add_log('Warning',iDict['filename'],infoCard, f"New {iDict['card_type']} VITEK card",'-')
    else:
        if upload:
            #print(f"Update [{upload}] -> Info")
            valLog.add_log('Info',iDict['filename'],infoCard, f"Update [{iDict['card_type']}] VITEK card",'-')
        else:
            #print(f"Update [{upload}] -> Warning")
            valLog.add_log('Warning',iDict['filename'],infoCard, f"Update [{iDict['card_type']}] VITEK card",'-')

    djVitekCard.orgbatch_id = OrgBatch
    djVitekCard.card_type = Dictionary.get(djVitekCard.Choice_Dictionary["card_type"],iDict['card_type'])
    if djVitekCard.card_type is None:
        valLog.add_log('Error',iDict['filename'],iDict['card_type'],'Vitek Card Type not correct','-')
        validStatus = False

    djVitekCard.card_code = iDict['card_code']
    djVitekCard.instrument = iDict['instrument']

    # in case cannot no date parsed from .pdf
    djVitekCard.expiry_date = iDict['expiry_date'] if 'expiry_date' in iDict.keys() else None
    djVitekCard.proc_date = iDict['processing_date'] if 'processing_date' in iDict.keys() else None
    djVitekCard.analysis_time = iDict['analysis_time'] if 'analysis_time' in iDict.keys() else None

    djVitekCard.clean_Fields()
    validDict = djVitekCard.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',iDict['filename'],k,validDict[k],'-')

    djVitekCard.VALID_STATUS = validStatus
    return(djVitekCard)

# ----------------------------------------------------------------------------------------------------
def imp_VitekID_fromDict(iDict,valLog,upload=False):
    """
    Create VITEK_ID instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    Barcode = VITEK_Card.get(iDict['card_barcode']) 
    if Barcode is None:
        if upload:
            validStatus = False
            valLog.add_log('Error',iDict['filename'],f"{iDict['card_barcode']} ({iDict['card_code']})",'[ID] VITEK card does not exists','-')
        else:
            valLog.add_log('Info',iDict['filename'],f"{iDict['card_barcode']} ({iDict['card_code']})",'New [ID] VITEK card will be uploaded','-')

    djVitekID = VITEK_ID.get(Barcode)
    if djVitekID is None:
        djVitekID = VITEK_ID()
        djVitekID.card_barcode = Barcode
        valLog.add_log('Info',iDict['filename'],f"{iDict['card_barcode']} ({iDict['card_code']})",'New VITEK ID','-')
    else:
        valLog.add_log('Info',iDict['filename'],f"{iDict['card_barcode']} ({iDict['card_code']})",'Update VITEK ID','-')

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
            if k == 'card_barcode':
                if upload:
                    valLog.add_log('Error',iDict['filename'],f"[ID] {k}",validDict[k],'-')
            else:
                valLog.add_log('Warning',iDict['filename'],f"[ID] {k}",validDict[k],'-')

    djVitekID.VALID_STATUS = validStatus
    return(djVitekID)


# ----------------------------------------------------------------------------------------------------
def imp_VitekAST_fromDict(iDict,valLog,upload=False):
    """
    Create VITEK_ID instance from a {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True
    Barcode = VITEK_Card.get(iDict['card_barcode']) 
    if Barcode is None:
        if upload:
            validStatus = False
            valLog.add_log('Error',iDict['filename'],f"{iDict['card_code']} ({iDict['card_barcode']})",'[AST] VITEK card does not exists','-')
        else:
            valLog.add_log('Info',iDict['filename'],f"{iDict['card_barcode']} ({iDict['card_code']})",'New [AST] VITEK card will be uploaded','-')

    DrugID = Drug.get(iDict['drug_name'])
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error',iDict['filename'],f"{iDict['drug_name']} ({iDict['card_barcode']})",'[AST] VITEK Drug does not Exists','-')

    if validStatus:
        djVitekAST = VITEK_AST.get(Barcode,DrugID,iDict['bp_source'],iDict['selected_organism'])
    else:
        djVitekAST = None
            
    if djVitekAST is None:
        djVitekAST = VITEK_AST()
        djVitekAST.card_barcode = Barcode
        djVitekAST.drug_id = DrugID
        djVitekAST.bp_source = iDict['bp_source']
        valLog.add_log('Info',iDict['filename'],f"{iDict['drug_name']} ({iDict['bp_source']})",'New [AST] VITEK','-')
    
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
            if k == 'card_barcode':
                if upload:
                    valLog.add_log('Error',iDict['filename'],f"[AST] {k}",validDict[k],'-')
            else:
                valLog.add_log('Warning',iDict['filename'],f"[AST] {k}",validDict[k],'-')

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
        valLog.add_log('Error','',f"{iDict['drug_name']} ",'Drug does not Exists','-')

    OrgBatchID = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatchID is None:
        validStatus = False
        valLog.add_log('Error','',f"{iDict['orgbatch_id']} ",'OrgBatchID does not Exists','-')

    RunID = Screen_Run.get(iDict['run_id'])
    if RunID is None:
        validStatus = False
        valLog.add_log('Error','',f"{iDict['run_id']} ",'RunID does not Exists','-')

    if validStatus:
        djMIC = MIC_COADD.get(OrgBatchID,DrugID,iDict['testplate_id'],iDict['testwell_id'])
    else:
        djMIC = None
            
    if djMIC is None:
        djMIC = MIC_COADD()
        djMIC.orgbatch_id = OrgBatchID
        djMIC.drug_id = DrugID
        djMIC.run_id = RunID
        djMIC.testplate_id = iDict['testplate_id']
        djMIC.testwell_id = iDict['testwell_id']
        valLog.add_log('Info','',f"{OrgBatchID} {DrugID} {iDict['testplate_id']}:{iDict['testwell_id']}",'New MIC ','-')
    
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
            valLog.add_log('Warning','',k,validDict[k],'-')

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
        valLog.add_log('Error','',f"{iDict['drug_name']} ",'Drug does not Exists','-')

    OrganismID = Organism.get(iDict['organism_id']) 
    if OrganismID is None:
        validStatus = False
        valLog.add_log('Error','',f"{iDict['organism_id']} ",'OrganismID does not Exists','-')

    if validStatus:
        djMIC = MIC_Pub.get(OrganismID,DrugID,iDict['source'])
    else:
        djMIC = None
            
    if djMIC is None:
        djMIC = MIC_Pub()
        djMIC.organism_id = OrganismID
        djMIC.drug_id = DrugID
        djMIC.source = iDict['source']
        valLog.add_log('Info','',f"{iDict['organism_id']} {iDict['drug_name']} {iDict['source']}",'New MIC ','-')
    
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
            valLog.add_log('Warning','',k,validDict[k],'-')

    djMIC.VALID_STATUS = validStatus
    
    return(djMIC)
