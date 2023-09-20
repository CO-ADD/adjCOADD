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
from dgene.utils.import_gene import imp_Sequence_fromDict, imp_FastQC_fromDict, imp_CheckM_fromDict,imp_IDSeq_fromDict
from dgene.utils.parse_wgs import split_BatchID_RunID, get_FastQC_Info, get_CheckM_Info, get_Kraken_Info, get_MLST_Info, get_GTDBTK_Info

from apputil.utils.data import listFolders
from apputil.utils import validation_log

# ----------------------------------------------------------------------------------------------------
def run_impProcess(lstDict,toDictFunc,procName,appuser=None,upload=False,OutputN=100):
# ----------------------------------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    nTotal = len(lstDict)
    nTime  = zData.Timer(nTotal)
    vLog = validation_log.Validation_Log(procName)
    logger.info(f"[{procName}] Processing {nTotal} for import")
    for f in lstDict:
        vLog.reset()
        djInst = toDictFunc(f,vLog)

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
    logger.info(f"[{procName}] [{nTotal-(nProc['Saved']+nProc['notValid'])}] {nTotal} {nProc}")

#-----------------------------------------------------------------------------------
def update_WGSCOADD_zAssembly(upload=False,uploaduser=None,OutputN=100):  
#-----------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    MicroOrgDB = "Z:/MicroOrgDB-Q5308"
    AssemblyBase = os.path.join(MicroOrgDB,"Sequence","WGS","02_Assembly")

    lstFastQC = []
    lstCheckM = []
    dictSequence = {}

    vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for dirHead in listFolders(AssemblyBase):
        zAssemblyFolder = os.path.join(AssemblyBase,dirHead)
        for BatchRunID in listFolders(zAssemblyFolder):
            dirAss = os.path.join(zAssemblyFolder,f"{BatchRunID}")
            if os.path.exists(dirAss):

                nProcessed = nProcessed + 1

                # Sequences -----------------------------
                seqDict={}
                seqDict['seq_name'] = BatchRunID
                seqDict['orgbatch_id'], seqDict['run_id'] = split_BatchID_RunID(seqDict['seq_name'])

                seqDict['seq_type'] = 'WGS'
                seqDict['seq_method'] = 'Illumina'
                seqDict['source'] = 'CO-ADD'
                seqDict['source_code'] = BatchRunID
                seqDict['source_link'] = f"RDM: {dirHead};{BatchRunID}"
                seqDict['reference'] = ""

                djSeq = imp_Sequence_fromDict(seqDict,vLog) 
                if djSeq.VALID_STATUS:
                    if upload:
                        djSeq.save(user=appuser)
                        seqDict['seq_id'] = djSeq
                else:
                    vLog.show(logTypes= ['Error'])

                dictSequence[BatchRunID] = seqDict

                if nProcessed%OutputN == 0:
                    print(f"[WGS-Assembly] {nProcessed:5d} {dirHead} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                sDict = {'seq_name':BatchRunID}

                # FastQC -----------------------------
                lFastQC=get_FastQC_Info(dirAss,seqDict['orgbatch_id'],seqDict['run_id'])

                for row in lFastQC:
                    djFQc = imp_FastQC_fromDict(row,vLog, objSeq = djSeq)
                    if djFQc.VALID_STATUS:
                        if upload:
                            djFQc.save(user=appuser)
                        else:
                            vLog.show(logTypes= ['Error'])

                    lstFastQC.append(dict(sDict,**row))

                # CheckM -----------------------------
                lCheckM=get_CheckM_Info(dirAss,seqDict['orgbatch_id'],seqDict['run_id'])
                for row in lCheckM:
                    djCheckM = imp_CheckM_fromDict(row,vLog, objSeq = djSeq)
                    #print(djCheckM.VALID_STATUS)
                    if djCheckM.VALID_STATUS:
                        
                        if upload:
                            djCheckM.save(user=appuser)
                        else:
                            vLog.show(logTypes= ['Error'])
                    lstCheckM.append(dict(sDict,**row))

#-----------------------------------------------------------------------------------
def update_WGSCOADD_FastA(upload=False,uploaduser=None,OutputN=100):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    MicroOrgDB = "Z:/MicroOrgDB-Q5308"
    FastABase = os.path.join(MicroOrgDB,"Sequence","WGS","03_FastA")

    lstFastQC = []
    lstCheckM = []
    dictSequence = {}

    vLog = validation_log.Validation_Log('WGS-FastA')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for dirHead in listFolders(FastABase):
        zFastAFolder = os.path.join(FastABase,dirHead)
        for BatchRunID in listFolders(zFastAFolder):
            dirFA = os.path.join(zFastAFolder,f"{BatchRunID}")
            if os.path.exists(dirFA):

                nProcessed = nProcessed + 1

                # Sequences -----------------------------
                seqDict={}
                seqDict['seq_name'] = BatchRunID
                seqDict['orgbatch_id'], seqDict['run_id'] = split_BatchID_RunID(seqDict['seq_name'])

                seqDict['seq_type'] = 'WGS'
                seqDict['seq_method'] = 'Illumina'
                seqDict['seq_file'] = 'Contigs'
                seqDict['source'] = 'CO-ADD'
                seqDict['source_code'] = BatchRunID
                seqDict['source_link'] = f"RDM: {dirHead};{BatchRunID}"
                seqDict['reference'] = ""

                djSeq = imp_Sequence_fromDict(seqDict,vLog) 
                if djSeq.VALID_STATUS:
                    if upload:
                        djSeq.save(user=appuser)
                        seqDict['seq_id'] = djSeq
                else:
                    vLog.show(logTypes= ['Error'])

                #dictSequence[BatchRunID] = seqDict

                seqDict['seq_id'] = djSeq
                if nProcessed%OutputN == 0:
                    print(f"[WGS-FastA] {nProcessed:5d} {dirHead} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                # Kraken2 -----------------------------
                lKraken=get_Kraken_Info(dirFA,seqDict['orgbatch_id'],seqDict['run_id'])
                #sDict = {'OrgBatch_ID':fa['orgbatchid'], 'Run_ID':fa['runid'], 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastA','Source':'CO-ADD'}
                idLst = []
                for v in lKraken:
                    # Formatted String for ArrayField
                    idLst.append(f"[{v['pct']:.1f} pct] {v['org_name']} ({v['tax_id']})")
                #print(idLst)
                seqDict['kraken_organisms'] = idLst

                # MLST -----------------------------
                lMLST=get_MLST_Info(dirFA,seqDict['orgbatch_id'],seqDict['run_id'])
                if len(lMLST) >0:
                    seqDict['mlst_scheme'] = lMLST[0]['mlst_scheme']
                    seqDict['mlst_seqtype'] = lMLST[0]['mlst_seqtype']
                    seqDict['mlst_alleles'] = lMLST[0]['mlst_alleles']

                # GTDBTK -----------------------------
                lGT=get_GTDBTK_Info(dirFA,seqDict['orgbatch_id'],seqDict['run_id'])
                if len(lGT) >0:
                    seqDict['gtdbtk_class'] = lGT[0]['gtdbtk_class']
                    seqDict['gtdbtk_fastani'] = f"{lGT[0]['gtdbtk_fastani_ref']} ({lGT[0]['gtdbtk_fastani_ani']})"

                #print(seqDict)
                djIDSeq = imp_IDSeq_fromDict(seqDict,vLog, objSeq = djSeq)
                if djIDSeq.VALID_STATUS:
                    #print(djIDSeq.VALID_STATUS)
                    if upload:
                        djIDSeq.save(user=appuser)
                    else:
                        vLog.show(logTypes= ['Error'])

