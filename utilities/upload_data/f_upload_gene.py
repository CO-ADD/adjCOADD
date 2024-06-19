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
#from dgene.utils.import_gene import (imp_Sequence_fromDict)
#from dgene.utils.parse_wgs import ()
from dgene.utils.upload_gene import (get_RDM, split_BatchID_RunID, get_subdir,
                                     upload_Trim, upload_CheckM, upload_FastA, upload_AMR)
from apputil.utils.data import listFolders
from apputil.utils import validation_log

#-----------------------------------------------------------------------------------
def update_WGSCOADD_Trim(upload=False,uploaduser=None,OutputN=100):  
#-----------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    TrimBase = os.path.join(RDM['base'],RDM['fastq_trim'])

    lstFastQC = []
    lstCheckM = []
    dictSequence = {}

    vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for subDir in listFolders(TrimBase):
        TrimFolder = os.path.join(TrimBase,subDir)
        for BatchRunID in listFolders(TrimFolder):
            dirTrim = os.path.join(TrimFolder,f"{BatchRunID}")
            if os.path.exists(dirTrim):

                nProcessed = nProcessed + 1
                OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)

                upload_Trim(OrgBatchID, RunID, dirTrim, vLog, upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-Assembly] {nProcessed:5d} {subDir} {OrgBatchID} {RunID} ")


                # nProcessed = nProcessed + 1

                # # Sequences -----------------------------
                # seqDict = update_GenomeSequence(subDir,BatchRunID,'WGS','Illumina','CO-ADD','Contigs','',vLog,upload=upload,uploaduser=uploaduser)

                # dictSequence[BatchRunID] = seqDict

                # if nProcessed%OutputN == 0:
                #     print(f"[WGS-Assembly] {nProcessed:5d} {subDir} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                # sDict = {'seq_name':BatchRunID}

                # # FastQC -----------------------------
                # lFastQC=get_FastQC_Info(dirTrim,seqDict['orgbatch_id'],seqDict['run_id'])

                # for row in lFastQC:
                #     djFQc = imp_FastQC_fromDict(row,vLog, objSeq = seqDict['seq_id'])
                #     if djFQc.VALID_STATUS:
                #         if upload:
                #             djFQc.save(user=appuser)
                #         else:
                #             vLog.show(logTypes= ['Error'])

                #     lstFastQC.append(dict(sDict,**row))


#-----------------------------------------------------------------------------------
def update_WGSCOADD_Assembly_single(OrgBatchID, RunID, upload=False,uploaduser=None):  
#-----------------------------------------------------------------------------

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    AssemblyBase = os.path.join(RDM['base'],RDM['assembly'])

    vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    subDir = get_subdir(OrgBatchID)
    dirAss = os.path.join(AssemblyBase,subDir,f"{OrgBatchID}_{RunID}")
    if os.path.exists(dirAss):
        upload_CheckM(OrgBatchID, RunID, dirAss, vLog, upload=upload,uploaduser=uploaduser)
        if upload:
            print(f"[UPLOADED]  {subDir} {OrgBatchID} {RunID} ")
        else:
            print(f"[PARSED  ]  {subDir} {OrgBatchID} {RunID} ")
    else:
        print(f"[FAILED  ]  WGS-Assembly: {dirAss} not found ")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_Assembly(upload=False,uploaduser=None,OutputN=100):  
#-----------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    AssemblyBase = os.path.join(RDM['base'],RDM['assembly'])

    vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for subDir in listFolders(AssemblyBase):
        zAssemblyFolder = os.path.join(AssemblyBase,subDir)
        for BatchRunID in listFolders(zAssemblyFolder):
            dirAss = os.path.join(zAssemblyFolder,f"{BatchRunID}")
            if os.path.exists(dirAss):

                nProcessed = nProcessed + 1
                OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)

                upload_CheckM(OrgBatchID, RunID, dirAss, vLog, upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-Assembly] {nProcessed:5d} {subDir} {OrgBatchID} {RunID} ")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_FastA_single(OrgBatchID, RunID, upload=False,uploaduser=None):  
#-----------------------------------------------------------------------------

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    FastABase = os.path.join(RDM['base'],RDM['fasta'])

    vLog = validation_log.Validation_Log('WGS-FastA')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    subDir = get_subdir(OrgBatchID)
    dirFA = os.path.join(FastABase,subDir,f"{OrgBatchID}_{RunID}")
    if os.path.exists(dirFA):
        upload_FastA(OrgBatchID, RunID, dirFA, vLog, upload=upload,uploaduser=uploaduser)
        if upload:
            print(f"[UPLOADED]  {subDir} {OrgBatchID} {RunID} ")
        else:
            print(f"[PARSED  ]  {subDir} {OrgBatchID} {RunID} ")
    else:
        print(f"[FAILED  ]  WGS-FastA: {dirFA} not found ")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_FastA(upload=False,uploaduser=None,OutputN=100):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    FastABase = os.path.join(RDM['base'],RDM['fasta'])

    vLog = validation_log.Validation_Log('WGS-FastA')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for subDir in listFolders(FastABase):
        zFastAFolder = os.path.join(FastABase,subDir)
        for BatchRunID in listFolders(zFastAFolder):
            dirFA = os.path.join(zFastAFolder,f"{BatchRunID}")
            if os.path.exists(dirFA):

                nProcessed = nProcessed + 1
                OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)

                upload_FastA(OrgBatchID, RunID, dirFA, vLog, upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-FastA] {nProcessed:5d} {subDir} {OrgBatchID} {RunID} ")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_AMR_single(OrgBatchID, RunID, Methods= ['AMR Finder'], upload=False,uploaduser=None):  
#-----------------------------------------------------------------------------

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    FastABase = os.path.join(RDM['base'],RDM['fasta'])

    vLog = validation_log.Validation_Log('WGS-FastA')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    subDir = get_subdir(OrgBatchID)
    dirFA = os.path.join(FastABase,subDir,f"{OrgBatchID}_{RunID}")
    if os.path.exists(dirFA):
        upload_AMR(OrgBatchID, RunID, dirFA, vLog, Methods, upload=upload,uploaduser=uploaduser)
        if upload:
            print(f"[UPLOADED]  {subDir} {OrgBatchID} {RunID} {Methods}")
        else:
            print(f"[PARSED  ]  {subDir} {OrgBatchID} {RunID} {Methods}")
    else:
        print(f"[FAILED  ]  WGS-AMR: {dirFA} not found ")


#-----------------------------------------------------------------------------------
def update_WGSCOADD_AMR(Methods= ['AMR Finder'], upload=False,uploaduser=None,OutputN=10):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    RDM = get_RDM("Z:/MicroOrgDB-Q5308")
    FastABase = os.path.join(RDM['base'],RDM['fasta'])

    vLog = validation_log.Validation_Log('WGS-FastA')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for subDir in listFolders(FastABase):
        zFastAFolder = os.path.join(FastABase,subDir)
        for BatchRunID in listFolders(zFastAFolder):
            dirFA = os.path.join(zFastAFolder,f"{BatchRunID}")
            if os.path.exists(dirFA):

                nProcessed = nProcessed + 1
                OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)

                upload_AMR(OrgBatchID, RunID, dirFA, vLog, Methods, upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-AMR] {nProcessed:5d} {subDir} {OrgBatchID} {RunID} {Methods} ")
