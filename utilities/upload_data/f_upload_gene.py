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
from dgene.utils.import_gene import (imp_Sequence_fromDict, 
                                     imp_FastQC_fromDict, imp_CheckM_fromDict,
                                     imp_IDSeq_fromDict, 
                                     imp_Gene_fromDict, imp_AMRGenotype_fromDict)
from dgene.utils.parse_wgs import (split_BatchID_RunID, 
                                   get_FastQC_Info, get_CheckM_Info, 
                                   get_Kraken_Info, get_MLST_Info, get_GTDBTK_Info, 
                                   get_AMRFinder_Info, get_Abricate_Info)

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
def update_GenomeSequence(SubDir,BatchRunID,SeqType,SeqMethod,SeqSource,SeqFile,SeqRef,vLog,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------
    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    seqDict={}
    seqDict['seq_name'] = BatchRunID
    seqDict['orgbatch_id'], seqDict['run_id'] = split_BatchID_RunID(seqDict['seq_name'])

    seqDict['seq_type'] = SeqType
    seqDict['seq_method'] = SeqMethod
    seqDict['source'] = SeqSource
    seqDict['source_code'] = BatchRunID
    seqDict['source_link'] = f"RDM: {SubDir};{BatchRunID}"
    seqDict['seq_file'] = SeqFile
    seqDict['reference'] = SeqRef

    djSeq = imp_Sequence_fromDict(seqDict,vLog) 
    if djSeq.VALID_STATUS:
        if upload:
            djSeq.save(user=appuser)    
    else:
        vLog.show(logTypes= ['Error'])
    seqDict['seq_id'] = djSeq
    return(seqDict)

#-----------------------------------------------------------------------------------
def update_Gene(GeneDict,vLog,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------
    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    djGene = imp_Gene_fromDict(GeneDict,vLog)

    if djGene.VALID_STATUS:
        if upload:
            djGene.save(user=appuser)
    else:
        vLog.show(logTypes= ['Error'])
    
    GeneDict['gene_id'] = djGene
    return(GeneDict)

#-----------------------------------------------------------------------------------
def update_WGSCOADD_Trim(upload=False,uploaduser=None,OutputN=100):  
#-----------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    MicroOrgDB = "Z:/MicroOrgDB-Q5308"
    TrimBase = os.path.join(MicroOrgDB,"Sequence","WGS","01_FastQ_Trim")

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

                # Sequences -----------------------------
                seqDict = update_GenomeSequence(subDir,BatchRunID,'WGS','Illumina','CO-ADD','Contigs','',vLog,upload=upload,uploaduser=uploaduser)

                dictSequence[BatchRunID] = seqDict

                if nProcessed%OutputN == 0:
                    print(f"[WGS-Assembly] {nProcessed:5d} {subDir} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                sDict = {'seq_name':BatchRunID}

                # FastQC -----------------------------
                lFastQC=get_FastQC_Info(dirTrim,seqDict['orgbatch_id'],seqDict['run_id'])

                for row in lFastQC:
                    djFQc = imp_FastQC_fromDict(row,vLog, objSeq = seqDict['seq_id'])
                    if djFQc.VALID_STATUS:
                        if upload:
                            djFQc.save(user=appuser)
                        else:
                            vLog.show(logTypes= ['Error'])

                    lstFastQC.append(dict(sDict,**row))


#-----------------------------------------------------------------------------------
def update_WGSCOADD_Assembly(upload=False,uploaduser=None,OutputN=100):  
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

    for subDir in listFolders(AssemblyBase):
        zAssemblyFolder = os.path.join(AssemblyBase,subDir)
        for BatchRunID in listFolders(zAssemblyFolder):
            dirAss = os.path.join(zAssemblyFolder,f"{BatchRunID}")
            if os.path.exists(dirAss):

                nProcessed = nProcessed + 1

                # Sequences -----------------------------
                seqDict = update_GenomeSequence(subDir,BatchRunID,'WGS','Illumina','CO-ADD','Contigs','',vLog,upload=upload,uploaduser=uploaduser)

                dictSequence[BatchRunID] = seqDict

                if nProcessed%OutputN == 0:
                    print(f"[WGS-Assembly] {nProcessed:5d} {subDir} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                sDict = {'seq_name':BatchRunID}

                # CheckM -----------------------------
                lCheckM=get_CheckM_Info(dirAss,seqDict['orgbatch_id'],seqDict['run_id'])
                for row in lCheckM:
                    djCheckM = imp_CheckM_fromDict(row,vLog, objSeq = seqDict['seq_id'])
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

    for subDir in listFolders(FastABase):
        zFastAFolder = os.path.join(FastABase,subDir)
        for BatchRunID in listFolders(zFastAFolder):
            dirFA = os.path.join(zFastAFolder,f"{BatchRunID}")
            if os.path.exists(dirFA):

                nProcessed = nProcessed + 1

                # Sequences -----------------------------
                seqDict = update_GenomeSequence(subDir,BatchRunID,'WGS','Illumina','CO-ADD','Contigs','',vLog,upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-FastA] {nProcessed:5d} {subDir} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

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
                djIDSeq = imp_IDSeq_fromDict(seqDict,vLog, objSeq = seqDict['seq_id'])
                if djIDSeq.VALID_STATUS:
                    #print(djIDSeq.VALID_STATUS)
                    if upload:
                        djIDSeq.save(user=appuser)
                    else:
                        vLog.show(logTypes= ['Error'])

                

#-----------------------------------------------------------------------------------
def update_WGSCOADD_AMR(Methods= ['AMR Finder'], upload=False,uploaduser=None,OutputN=10):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProcessed = 0

    MicroOrgDB = "Z:/MicroOrgDB-Q5308"
    FastABase = os.path.join(MicroOrgDB,"Sequence","WGS","03_FastA")

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

                # Sequences -----------------------------
                seqDict = update_GenomeSequence(subDir,BatchRunID,'WGS','Illumina','CO-ADD','Contigs','',vLog,upload=upload,uploaduser=uploaduser)

                if nProcessed%OutputN == 0:
                    print(f"[WGS-AMR] {nProcessed:5d} {subDir} {seqDict['orgbatch_id']} {seqDict['run_id']} ")

                # AMR Finder -----------------------------
                if 'AMR Finder' in Methods:
                    lAmrFinder=get_AMRFinder_Info(dirFA,seqDict['orgbatch_id'],seqDict['run_id'])
                    for row in lAmrFinder:

                        #row['source'] = "AMR Finder"
                        djGene = imp_Gene_fromDict(row,vLog)
                        row['gene_id'] = djGene
                        if djGene.VALID_STATUS:
                            if upload:
                                djGene.save(user=appuser)
                                row['gene_id'] = djGene
                        else:
                            vLog.show(logTypes= ['Error'])

                        row['seq_id'] = seqDict['seq_id']
                        djAMRgt = imp_AMRGenotype_fromDict(row,vLog)
                        if djAMRgt.VALID_STATUS:
                            if upload:
                                djAMRgt.save(user=appuser)
                        else:
                            vLog.show(logTypes= ['Error'])

                # Abricate CARD -----------------------------
                if 'Abricate card' in Methods:
                    lAbCard=get_Abricate_Info(dirFA,seqDict['orgbatch_id'],seqDict['run_id'],DB='card')
                    for row in lAbCard:

                        print(row)
                        djGene = imp_Gene_fromDict(row,vLog)
                        #print(f"[{len(djGene.gene_name):3d}] {djGene.gene_name}")
                        row['gene_id'] = djGene
                        if djGene.VALID_STATUS:
                            if upload:
                                djGene.save(user=appuser)
                                row['gene_id'] = djGene
                        else:
                            vLog.show(logTypes= ['Error'])
                        print("...")
                        row['seq_id'] = seqDict['seq_id']
                        djAMRgt = imp_AMRGenotype_fromDict(row,vLog)
                        if djAMRgt.VALID_STATUS:
                            if upload:
                                djAMRgt.save(user=appuser)
                        else:
                            vLog.show(logTypes= ['Error'])
