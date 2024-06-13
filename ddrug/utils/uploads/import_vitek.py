import os
from pathlib import Path
# from django_rdkit.models import *
# from django_rdkit.config import config
# from django.conf import settings

from dorganism.models import Organism_Batch
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST
from dscreen.models import Screen_Run
# from ddrug.utils.molecules import *

from apputil.models import ApplicationUser, Dictionary
# from apputil.utils.data import *

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
 
