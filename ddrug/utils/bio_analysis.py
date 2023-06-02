import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

#import django

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from dorganism.utils.utils import reformat_OrganismID, reformat_OrgBatchID
from ddrug.models import Drug, MIC_COADD, MIC_Pub, Breakpoint
from ddrug.utils.bio_data import split_XC50


#-------------------------------------------------------------------------------------------
def get_BreakPoint_Object(DrugName,OrgName,BPType):
    """
    Get the Breakpoint object based on Organism Name and BP Type (MIC/Zone)
    """
#-------------------------------------------------------------------------------------------

    objDrug = Drug.get(DrugName)
    objTax  = Taxonomy.get(OrgName)

    bpSource = ['EUCAST','CLSI']

    if objDrug and objTax:
        lstBP = Breakpoint.objects.filter(drug_id=objDrug,bp_type=BPType)
        TaxLineage = objTax.lineage
        TaxNameLst = OrgName.split(' ')

        selBP = (0,None)
        notBP = (0,None)
        for objBP in lstBP:
            _notRank = str(objBP.notorg_rank)
            _orgRank = str(objBP.org_rank)
            if objBP.bp_source == 'EUCAST':
                if _orgRank == 'Specie' and objBP.org_name == OrgName and selBP[0] < 10:
                    selBP = (10,objBP)
                elif _orgRank == 'Genus' and objBP.org_name in TaxNameLst and selBP[0] < 9:
                    selBP = (9,objBP)
                elif _orgRank == 'Family' and objBP.org_name in TaxLineage and selBP[0] < 8:
                    selBP = (8,objBP)

                # if _notRank == 'Specie' and objBP.notorg_name != OrgName and notBP[0] < 10:
                #     notBP = (10,objBP)
                # elif _notRank == 'Genus' and objBP.notorg_name not in TaxNameLst and notBP[0] < 9:
                #     notBP = (9,objBP)
                # elif _notRank == 'Family' and objBP.notorg_name not in TaxLineage and notBP[0] < 8:
                #     notBP = (8,objBP)

            if objBP.bp_source == 'CLSI':
                if _orgRank == 'Specie' and objBP.org_name == OrgName and selBP[0] < 5:
                    selBP = (5,objBP)
                elif _orgRank == 'Genus' and objBP.org_name in TaxNameLst and selBP[0] < 4:
                    selBP = (4,objBP)
                elif _orgRank == 'Family' and objBP.org_name in TaxLineage and selBP[0] < 3:
                    selBP = (3,objBP)

                # if _notRank == 'Specie' and objBP.notorg_name != OrgName and notBP[0] < 5:
                #     notBP = (5,objBP)
                # elif _notRank == 'Genus' and objBP.notorg_name not in TaxNameLst and notBP[0] < 4:
                #     notBP = (4,objBP)
                # elif _notRank == 'Family' and objBP.notorg_name not in TaxLineage and notBP[0] < 3:
                #     notBP = (3,objBP)
        if selBP[0] > 0:
            return(selBP[1])
        else:
            if notBP[0] > 0:
                return(notBP[1])
            else:
                logger.info(f"[get_BreakPoint_Object] Unable to find BP for {DrugName} {OrgName} {BPType} from {len(lstBP)} entries")
                return(None)
    else:
        logger.info(f"[DEBUG] No Breakpoints found for {DrugName} [{objDrug}] against {OrgName} [{objTax}]")
    return(None)

#-----------------------------------------------------------------------------
def calc_bp(valAct,objBP):
    """
    Calculate the Breakpoint profile from MIC/Zone value and Breakpoint object
    """
#-----------------------------------------------------------------------------
    if objBP:
        _bpType = str(objBP.bp_type)
        _prefix,_val,_sval = split_XC50(valAct)
        if _bpType == 'MIC':
            _bp = 'I'
            if _val <= objBP.bp_sens_le :
                _bp = 'S'
            elif _val > objBP.bp_res_gt :
                _bp = 'R'
            return _bp
        if _bpType == 'Zone':
            _bp = 'I'
            if _val <= objBP.bp_sens_le :
                _bp = 'S'
            elif _val > objBP.bp_res_gt :
                _bp = 'R'
            return _bp
    return(None)

#-----------------------------------------------------------------------------
def calc_Breakpoint(DrugName,OrgName,BPType,DRVal):
    """
    Calculate the Breakpoint based on Drug Name, Organism Name, BP Type (MIC/Zone) and DRValue
    """
#-----------------------------------------------------------------------------
    djBP = get_BreakPoint_Object(DrugName,OrgName,BPType)
    if djBP:
        bp = calc_bp(DRVal,djBP)
        return(bp,djBP.info())
    return("-","no BP")
