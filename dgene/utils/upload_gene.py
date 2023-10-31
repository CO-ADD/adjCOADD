import os
from datetime import datetime
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

#from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
#from dscreen.models import Screen_Run
from dgene.models import Genome_Sequence,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM, Gene, AMR_Genotype
from dgene.utils.import_gene import (imp_Sequence_fromDict, 
                                     imp_FastQC_fromDict, imp_CheckM_fromDict,
                                     imp_IDSeq_fromDict, 
                                     imp_Gene_fromDict, imp_AMRGenotype_fromDict)
from dgene.utils.parse_wgs import (get_FastQC_Info, get_CheckM_Info, 
                                   get_Kraken_Info, get_MLST_Info, get_GTDBTK_Info, 
                                   get_AMRFinder_Info, get_Abricate_Info)
 
from apputil.models import ApplicationUser, Dictionary
#from apputil.utils.data import *

#-----------------------------------------------------------------------------
def get_RDM(MicroOrgDB):
#-----------------------------------------------------------------------------
    RDM = {
    'rdm' : MicroOrgDB,
    'base' : os.path.join(MicroOrgDB,"Sequence","WGS"),
    'fastq':'01_FastQ',
    'fastq_trim':'01_FastQ_Trim',
    'assembly':'02_Assembly',
    'fasta':'03_FastA',
    }
    return(RDM)

#-----------------------------------------------------------------------------
def split_BatchID_RunID(batch_run_id):
#-----------------------------------------------------------------------------
    arrStr = batch_run_id.split("_")
    batchID = '_'.join(arrStr[0:3])
    runID = '_'.join(arrStr[3:])
    return batchID, runID

#-----------------------------------------------------------------------------
def get_subdir(OrgBatchID,binsize=200):
#-----------------------------------------------------------------------------
    """
     Gets the SubFolder name based on the XX_NNNN with splits into 200
      GN_0000, GN_0200, GN_0400, ... ,GN_1200, GN_1400, GN_1600
    """
    _org = OrgBatchID.split('_')
    return(f"{_org[0]}_{int(int(_org[1])/binsize)*binsize:04d}")

#-----------------------------------------------------------------------------------
def gen_SeqDict(OrgBatchID,RunID,SeqType,SeqMethod,Source):
#-----------------------------------------------------------------------------------
    SeqDict = {
        'seq_name'   : f"{OrgBatchID}_{RunID}",
        'orgbatch_id': OrgBatchID,
        'run_id'     : RunID,
        'seq_type'   : SeqType,
        'seq_method' : SeqMethod,
        'source'     : Source,
        'source_code': f"{OrgBatchID}_{RunID}",
        'source_link': f"RDM {SeqType}: {OrgBatchID}_{RunID}",
        'seq_file'   : 'Contigs',
        'reference'  : ''
    }
    return(SeqDict)

#-----------------------------------------------------------------------------------
def upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
                          vLog,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------
    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    djSeq = imp_Sequence_fromDict(SeqDict,vLog) 
    if djSeq.VALID_STATUS:
        if upload:
            djSeq.save(user=appuser)    
    else:
        vLog.show(logTypes= ['Error'])
    return(djSeq)

#-----------------------------------------------------------------------------------
def upload_Gene(GeneDict,vLog,upload=False,uploaduser=None):
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
    
    return(djGene)

#-----------------------------------------------------------------------------------
def upload_Trim(OrgBatchID, RunID, TrimDir, vLog, upload=False,uploaduser=None,verbose=False):  
#-----------------------------------------------------------------------------
    lstFastQC = []
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(TrimDir):

        SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')

        # Sequences -----------------------------
        SeqDict['seq_id'] = upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
                                        vLog,upload=upload,uploaduser=uploaduser)

        if verbose:
            print(f"[WGS-Trim] {TrimDir} {OrgBatchID} {RunID} ")

        sDict = {'seq_name':f"{OrgBatchID}_{RunID}"}

        # FastQC -----------------------------
        lFastQC=get_FastQC_Info(TrimDir,OrgBatchID, RunID,)

        for row in lFastQC:
            djFQc = imp_FastQC_fromDict(row,vLog, objSeq = SeqDict['seq_id'])
            if djFQc.VALID_STATUS:
                if upload:
                    djFQc.save(user=appuser)
                else:
                    vLog.show(logTypes= ['Error'])

            lstFastQC.append(dict(sDict,**row))

#-----------------------------------------------------------------------------------
def upload_CheckM(OrgBatchID, RunID, AssemblyDir, vLog, upload=False,uploaduser=None,verbose=False):  
#-----------------------------------------------------------------------------

    lstCheckM = []

    #vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(AssemblyDir):

        SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')

        # Sequences -----------------------------
        SeqDict['seq_id'] = upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
                                        vLog,upload=upload,uploaduser=uploaduser)

        if verbose:
            print(f"[WGS-Assembly] {AssemblyDir} {OrgBatchID} {RunID} ")

        sDict = {'seq_name':f"{OrgBatchID}_{RunID}"}

        # CheckM -----------------------------
        lCheckM=get_CheckM_Info(AssemblyDir,SeqDict['orgbatch_id'],SeqDict['run_id'])
        for row in lCheckM:
            djCheckM = imp_CheckM_fromDict(row,vLog, objSeq = SeqDict['seq_id'])
            #print(djCheckM.VALID_STATUS)
            if djCheckM.VALID_STATUS:
                
                if upload:
                    djCheckM.save(user=appuser)
                else:
                    vLog.show(logTypes= ['Error'])
            lstCheckM.append(dict(sDict,**row))
    return(lstCheckM)    

#-----------------------------------------------------------------------------------
def upload_FastA(OrgBatchID, RunID, FastaDir, vLog, upload=False,uploaduser=None,verbose=False):  
#-----------------------------------------------------------------------------

    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(FastaDir):

        SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')
        # Sequences -----------------------------
        SeqDict['seq_id'] = upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
                                        vLog,upload=upload,uploaduser=uploaduser)

        if verbose:
            print(f"[WGS-Fasta] {FastaDir} {OrgBatchID} {RunID} ")

        #sDict = {'seq_name':f"{OrgBatchID}_{RunID}"}

        # Kraken2 -----------------------------
        lKraken=get_Kraken_Info(FastaDir,OrgBatchID, RunID)
        #sDict = {'OrgBatch_ID':fa['orgbatchid'], 'Run_ID':fa['runid'], 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastA','Source':'CO-ADD'}
        idLst = []
        for v in lKraken:
            # Formatted String for ArrayField
            idLst.append(f"{v['org_name']} ({v['tax_id']}) [{v['pct']:.1f} pct]")
        #print(idLst)
        SeqDict['kraken_organisms'] = idLst

        # MLST -----------------------------
        lMLST=get_MLST_Info(FastaDir,OrgBatchID, RunID)
        if len(lMLST) >0:
            SeqDict['mlst_scheme'] = lMLST[0]['mlst_scheme']
            SeqDict['mlst_seqtype'] = lMLST[0]['mlst_seqtype']
            SeqDict['mlst_alleles'] = lMLST[0]['mlst_alleles']

        # GTDBTK -----------------------------
        lGT=get_GTDBTK_Info(FastaDir,OrgBatchID, RunID)
        if len(lGT) >0:
            SeqDict['gtdbtk_class'] = lGT[0]['gtdbtk_class']
            SeqDict['gtdbtk_fastani'] = f"{lGT[0]['gtdbtk_fastani_ref']} ({lGT[0]['gtdbtk_fastani_ani']})"

        #print(seqDict)
        djIDSeq = imp_IDSeq_fromDict(SeqDict, vLog, objSeq = SeqDict['seq_id'])
        if djIDSeq.VALID_STATUS:
            #print(djIDSeq.VALID_STATUS)
            if upload:
                djIDSeq.save(user=appuser)
            else:
                vLog.show(logTypes= ['Error'])

#-----------------------------------------------------------------------------------
def upload_AMR(OrgBatchID, RunID, FastaDir, vLog, Methods= ['AMR Finder'], upload=False,uploaduser=None,verbose=False):
#-----------------------------------------------------------------------------------

    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(FastaDir):

        SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')
        # Sequences -----------------------------
        SeqDict['seq_id'] = upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
                                        vLog,upload=upload,uploaduser=uploaduser)

        if verbose:
            print(f"[WGS-AMR] {FastaDir} {OrgBatchID} {RunID} ")
                    # Sequences -----------------------------

        # AMR Finder -----------------------------
        if 'AMR Finder' in Methods:
            lAmrFinder=get_AMRFinder_Info(FastaDir,OrgBatchID, RunID)
            for row in lAmrFinder:

                row['gene_id'] = upload_Gene(row,vLog,upload=upload,uploaduser=uploaduser)
                row['seq_id'] = SeqDict['seq_id']

                djAMRgt = imp_AMRGenotype_fromDict(row,vLog)
                if djAMRgt.VALID_STATUS:
                    if upload:
                        djAMRgt.save(user=appuser)
                else:
                    vLog.show(logTypes= ['Error'])

        # Abricate CARD -----------------------------
        if 'Abricate card' in Methods:
            lAbCard=get_Abricate_Info(FastaDir,OrgBatchID, RunID,DB='card')
            for row in lAbCard:

                row['gene_id'] = upload_Gene(row,vLog,upload=upload,uploaduser=uploaduser)
                row['seq_id'] = SeqDict['seq_id']

                djAMRgt = imp_AMRGenotype_fromDict(row,vLog)
                if djAMRgt.VALID_STATUS:
                    if upload:
                        djAMRgt.save(user=appuser)
                else:
                    vLog.show(logTypes= ['Error'])