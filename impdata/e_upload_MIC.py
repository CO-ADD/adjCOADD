import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

import django
from oraCastDB import oraCastDB
from zUtils import zData

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, MIC_COADD, MIC_Pub
from apputil.utils import Validation_Log

#-----------------------------------------------------------------------------------
def reformat_OrganismID(OrgID):
#-----------------------------------------------------------------------------------
    xStr = OrgID.split("_")
    return(f"{xStr[0]}_{int(xStr[1]):04d}")

#-----------------------------------------------------------------------------------
def update_MICPub_ora(upload=False,uploaduser=None,OutputN=100):
#-----------------------------------------------------------------------------------
    micSQL = """
    Select Drug_ID, Drug_Code, Drug_Name, Organism_ID, Organism_Name,
        MIC, MIC_High, MIC_List, MIC_Low, MIC_Median, MIC_Mode, nMIC,
        MIC_Prefix, MIC_Value, MIC_Unit, 
        BP_Profile, BP_Res_GT, BP_Sen_LE, BP_Source,
        Media, Plate_Material, Plate_Size,
        Source, Source_ID, Source_Type
    From MIC
        Where Source <> 'CO-ADD'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[MIC-Pub] ... ")
    micLst = OrgDB.get_dict_list(micSQL)
    nTotal = len(micLst)
    logger.info(f"[MIC-Pub] {nTotal} ")
    OrgDB.close()

    vLog = Validation_Log('MIC-Pub')
    nTime  = zData.Timer(nTotal)
    nProcessed = 0

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0

    for mic in micLst:
        mic['ORGANISM_ID'] = reformat_OrganismID(mic['ORGANISM_ID'])
        djMIC = MIC_Pub.check_from_dict(mic,vLog)
        djMIC.clean_Fields()
        validDict = djMIC.validate()

        if validDict:
            logger.info(f"{mic['ORGANISM_ID']} {mic['DRUG_NAME']} {mic['SOURCE']} {validDict} ")

        # --- Upload ---------------------------------------------------------
        nProcessed = nProcessed + 1
        if djMIC.VALID_STATUS:
            if upload:
                logger.info(f" -> {djMIC} as {appuser}")
                djMIC.save(user=appuser)
                nProc['Saved'] = nProc['Saved'] + 1
            else:
                if nProcessed%OutputN == 0:
                    eTime,sTime = nTime.remains(nProcessed)
                    logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >r {djMIC} ")
        else:
            nProc['notValid'] = nProc['notValid'] + 1
    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[MIC-Pub] [{nTotal-(nProc['Saved']+nProc['notValid'])}] {nTotal} {nProc}")

