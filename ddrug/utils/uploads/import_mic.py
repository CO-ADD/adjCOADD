import os
from pathlib import Path
# from django_rdkit.models import *
# from django_rdkit.config import config
# from django.conf import settings

from dorganism.models import Organism, Organism_Batch
from ddrug.models import Drug, MIC_COADD, MIC_Pub, Breakpoint
from dscreen.models import Screen_Run
# from ddrug.utils.molecules import *

from apputil.models import ApplicationUser, Dictionary
# from apputil.utils.data import *

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