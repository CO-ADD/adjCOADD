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
from dorganism.utils.utils import reformat_OrganismID, reformat_OrgBatchID
from dgene.models import Gene,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM
from dgene.utils.import_gene import imp_FastQC_fromDict, imp_CheckM_fromDict,imp_IDSeq_fromDict
from dgene.utils.parse_wgs import parse_WGS_COADD
from apputil.utils import validation_log


def run_impProcess(lstDict,toDictFunc,procName,appuser=None,upload=False,OutputN=100):
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    nTotal = len(lstDict)
    nTime  = zData.Timer(nTotal)
    vLog = validation_log.Validation_Log('WGS-COADD')
    for f in lstDict:
        vLog.reset()
        djInst = toDictFunc(f,vLog)

        # --- Upload ---------------------------------------------------------
        nProcessed = nProcessed + 1
        if djInst.VALID_STATUS:
            if upload:
                if nProcessed%OutputN == 0:
                    eTime,sTime = nTime.remains(nProcessed)
                    logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} -Save-> {djInst} ")
                djInst.save(user=appuser)
                nProc['Saved'] = nProc['Saved'] + 1
            else:
                if nProcessed%OutputN == 0:
                    eTime,sTime = nTime.remains(nProcessed)
                    logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >Check< {djInst} ")
        else:
            vLog.show(logTypes= ['Error'])
            nProc['notValid'] = nProc['notValid'] + 1
    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[WGS-FastQC] [{nTotal-(nProc['Saved']+nProc['notValid'])}] {nTotal} {nProc}")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_zAssembly(upload=False,uploaduser=None,OutputN=100):
#-----------------------------------------------------------------------------------
    MicroOrgDB = "Z:/MicroOrgDB-Q5308"
    zAssemblyBase = os.path.join(MicroOrgDB,"Sequence","WGS","02_Assembly")
    lstFastQC, lstCheckM, lstKraken = parse_WGS_COADD(zAssemblyBase)

    vLog = validation_log.Validation_Log('WGS-COADD')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    # - FastQC
    # --------------------------------------------------------------------------------
    run_impProcess(lstFastQC,imp_FastQC_fromDict,"WGS-FastQC",appuser,upload=upload,OutputN=100)
    run_impProcess(lstCheckM,imp_CheckM_fromDict,"WGS-CheckM",appuser,upload=upload,OutputN=100)
    run_impProcess(lstKraken,imp_IDSeq_fromDict,"WGS-Kraken",appuser,upload=upload,OutputN=100)


