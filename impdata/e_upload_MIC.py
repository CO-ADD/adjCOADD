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
def reformat_OrgBatchID(OrgBatchID):
#-----------------------------------------------------------------------------------
    xStr = OrgBatchID.split(":")
    #print(xStr)
    return(f"{reformat_OrganismID(xStr[0])}_{int(xStr[1]):02d}")

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


#-----------------------------------------------------------------------------------
def update_MICCOADD_ora(RunID,upload=False,uploaduser=None,OutputN=100):
#-----------------------------------------------------------------------------------
    micSQL = """
    Select TestPlate_ID, TestWell_ID, 
      Compound_ID, Compound_Name, Compound_Code,
      Compound2_ID, Compound2_Name, Compound2_Code,
      Test_Strain, Plate_Size, Material, Test_Dye, Test_Additive, Test_Date,
      -- Media,
      DR, DR_Unit, DR2_Value, DR2_Unit, DMax, Data_Quality, Run_ID 
    From vDoseResponse
      Where Result_Type = 'MIC' 
    """
    micSQL += f" And Run_ID = '{RunID}'"
    CastDB = oraCastDB.openCastDB()
    logger.info(f"[MIC-COADD] ... ")
    micLst = CastDB.get_dict_list(micSQL)
    nTotal = len(micLst)
    logger.info(f"[MIC-COADD] {nTotal} ")
    CastDB.close()

    if nTotal>0:
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
            mic['ORGBATCH_ID'] = reformat_OrgBatchID(mic['TEST_STRAIN'])
            if mic['COMPOUND2_NAME']:
                mic['DRUG_NAME'] = mic['COMPOUND_NAME'] +'|'+ mic['COMPOUND2_NAME']
                if mic['DR2_VALUE']:
                    mic['MIC'] = mic['DR'] +'|'+ str(mic['DR2_VALUE'])
                else:
                    mic['MIC'] = mic['DR'] 
                if mic['DR2_UNIT']:
                    mic['MIC_UNIT'] = mic['DR_UNIT'] +'|'+ mic['DR2_UNIT']
                else:
                    mic['MIC_UNIT'] = mic['DR_UNIT']
            else:
                mic['DRUG_NAME'] = mic['COMPOUND_NAME']
                mic['MIC'] = mic['DR']
                mic['MIC_UNIT'] = mic['DR_UNIT']
            
            mic['MEDIA'] = None
            mic['PLATE_SIZE'] = mic['PLATE_SIZE'].replace('w','')
            mic['PLATE_MATERIAL'] = mic['MATERIAL']
            mic['DYE'] = mic['TEST_DYE']
            mic['ADDITIVE'] = mic['TEST_ADDITIVE']

            djMIC = MIC_COADD.check_from_dict(mic,vLog)
            djMIC.clean_Fields()
            validDict = djMIC.validate()

            if validDict:
                logger.info(f"{mic['ORGBATCH_ID']} {mic['DRUG_NAME']} {mic['RUN_ID']} {validDict} ")

            # --- Upload ---------------------------------------------------------
            nProcessed = nProcessed + 1
            if djMIC.VALID_STATUS:
                if upload:
                    if nProcessed%OutputN == 0:
                        eTime,sTime = nTime.remains(nProcessed)
                        logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} -> {djMIC} ")
                    djMIC.save(user=appuser)
                    nProc['Saved'] = nProc['Saved'] + 1
                else:
                    if nProcessed%OutputN == 0:
                        eTime,sTime = nTime.remains(nProcessed)
                        logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >r {djMIC} ")
            else:
                nProc['notValid'] = nProc['notValid'] + 1
        eTime,sTime = nTime.remains(nProcessed)
        logger.info(f"[MIC-COADD] [{nTotal-(nProc['Saved']+nProc['notValid'])}] {nTotal} {nProc}")
    else:
        logger.info(f"[MIC-COADD] [0 No records found]")
