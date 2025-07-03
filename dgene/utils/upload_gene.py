import os
import re
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
                                   get_AMRFinder_Info, get_Abricate_Info, get_RGI_Info)
 
from apputil.models import ApplicationUser, Dictionary
#from apputil.utils.data import *

#-----------------------------------------------------------------------------
class WGS_RDM():
#-----------------------------------------------------------------------------

    WGS_RDM_FOLDERS = {
        'fastq':'01_FastQ',
        'trim':'01_FastQ_Trim',
        'assembly':'02_Assembly',
        'fasta':'03_FastA',
    }

    def __init__(self,BaseWGS,OrgbatchID,RunID,SeqMethod='Illumina',Source='CO-ADD',valLog=None):
        

        self.wgs_base = BaseWGS
        self.seq_name = f"{OrgbatchID}_{RunID}"

        if '-' in OrgbatchID:
            self.orgbatch_id = OrgbatchID.split('-')[0]
        else:
            self.orgbatch_id = OrgbatchID
        self.run_id = RunID
        self.seq_type = 'WGS'
        self.seq_method =SeqMethod
        self.seq_dir = self.get_subdir(self.orgbatch_id)
        self.seq_file = 'Contigs'
        self.val_log = valLog
        self.kraken_organism =None

        self.seq_dict =  {
                    'seq_name'   : self.seq_name,
                    'orgbatch_id': self.orgbatch_id,
                    'run_id'     : self.run_id,
                    'seq_type'   : self.seq_type,
                    'seq_method' : self.seq_method,
                    'source'     : Source,
                    'source_code': self.seq_name,
                    'source_link': f"RDM {self.seq_type}: {self.orgbatch_id}_{self.run_id}",
                    'seq_file'   : self.seq_file,
                    'reference'  : ''
                }
       
        
        self.fastq_dir = os.path.join(self.wgs_base,self.WGS_RDM_FOLDERS['fastq'],self.seq_dir)
        self.trim_dir = os.path.join(self.wgs_base,self.WGS_RDM_FOLDERS['trim'],self.seq_dir,self.seq_name)
        self.assembly_dir = os.path.join(self.wgs_base,self.WGS_RDM_FOLDERS['assembly'],self.seq_dir,self.seq_name)
        self.fasta_dir = os.path.join(self.wgs_base,self.WGS_RDM_FOLDERS['fasta'],self.seq_dir,self.seq_name)

        # FastQ Files ----------------------
        self.fastq_files = {}
        self.trim_files = {}
        self.n_fastq = 0
        self.get_fastq_files()

        # FastQ Files ----------------------
        self.seq_id = Genome_Sequence.get(None,self.seq_name)


    #-------------------------------------------
    def __str__(self):
        return(self.seq_name)

    #-------------------------------------------
    def get_fastq_files(self):
        CHECK_FASTQ = {'pe1':'R1.fastq.gz','pe2':'R2.fastq.gz','se':'S.fastq.gz'}
        for ft in CHECK_FASTQ:
            #print(os.path.join(self.fastq_dir, f"{self.seq_name}_{CHECK_FASTQ[ft]}"))
            if os.path.isfile(os.path.join(self.fastq_dir, f"{self.seq_name}_{CHECK_FASTQ[ft]}")):
                self.fastq_files[ft] = f"{self.seq_name}_{CHECK_FASTQ[ft]}"

        CHECK_TRIM = {'pe1':'P1.fastq.gz','pe2':'P2.fastq.gz','se':'S.fastq.gz'}
        for ft in CHECK_TRIM:
            #print(os.path.join(self.fastq_dir, f"{self.seq_name}_{CHECK_FASTQ[ft]}"))
            if os.path.isfile(os.path.join(self.trim_dir, f"{self.seq_name}_{CHECK_TRIM[ft]}")):
                self.trim_files[ft] = f"{self.seq_name}_{CHECK_TRIM[ft]}"

        self.n_fastq = len(self.trim_files)+len(self.fastq_files)

    #-------------------------------------------
    @staticmethod
    def listFolders(Path):
        return([ name for name in os.listdir(Path) if os.path.isdir(os.path.join(Path, name)) ])

    #-------------------------------------------
    @staticmethod
    def get_subdir(OrgBatchID,binsize=200):
        """
        Gets the SubFolder name based on the XX_NNNN with splits into 200
        GN_0000, GN_0200, GN_0400, ... ,GN_1200, GN_1400, GN_1600
        """
        _org = OrgBatchID.split('_')
        return(f"{_org[0]}_{int(int(_org[1])/binsize)*binsize:04d}")

    #-----------------------------------------------------------------------------
    @staticmethod
    def split_BatchID_RunID(batch_run_id):
        arrStr = batch_run_id.split("_")
        batchID = '_'.join(arrStr[0:3])
        runID = '_'.join(arrStr[3:])
        return batchID, runID



    def upload_GenomeSequence(self,upload=False,uploaduser=None):
    #-----------------------------------------------------------------------------------
        # check user
        self.get_fastq_files()
        #print(self)
        if self.n_fastq > 0:
            appuser = None
            if uploaduser:
                appuser = ApplicationUser.get(uploaduser)

            self.seq_id = imp_Sequence_fromDict(self.seq_dict,self.val_log) 
            #print(f" [upload_GenomeSequence] {self.seq_id}")
            if self.seq_id:
                if self.seq_id.VALID_STATUS:
                    if upload:
                        self.seq_id.save(user=appuser)
            else:
                self.val_log.show(logTypes= ['Error'])
        else:
            print(f" [upload_GenomeSequence] No FastQ files found")

    #-----------------------------------------------------------------------------------
    def upload_CheckM(self, upload=False,uploaduser=None):
    #    vLog, upload=False,uploaduser=None,verbose=False):  
    #-----------------------------------------------------------------------------

        lstCheckM = []

        #vLog = validation_log.Validation_Log('WGS-Assembly')

        # check user
        appuser = None
        if uploaduser:
            appuser = ApplicationUser.get(uploaduser)

        if os.path.exists(self.assembly_dir):

            lCheckM=get_CheckM_Info(self.assembly_dir,self.orgbatch_id,self.run_id, 
                                    Assemblies = ['spades','shovill'], 
                                    outType = 'contigs_filtered', 
                                    Contamination_cutOff = 5.0)
            
            for row in lCheckM:
                if self.seq_id:
                    djCheckM = imp_CheckM_fromDict(row, self.val_log, self.seq_id)
                    # #print(djCheckM.VALID_STATUS)
                    if djCheckM.VALID_STATUS:
                        
                        if upload:
                            djCheckM.save(user=appuser)
                        else:
                            self.val_log.show(logTypes= ['Error'])
                #lstCheckM.append(dict(sDict,**row))

    #-----------------------------------------------------------------------------------
    def upload_FastA_ID(self, upload=False,uploaduser=None,verbose=False):  
    #-----------------------------------------------------------------------------

        appuser = None
        if uploaduser:
            appuser = ApplicationUser.get(uploaduser)

        if os.path.exists(self.fasta_dir):

            # Kraken2 -----------------------------
            lKraken=get_Kraken_Info(self.fasta_dir,self.orgbatch_id, self.run_id)
            self.seq_dict['kraken_organisms'] = []
            if len(lKraken) >0 :
                self.kraken_organism = lKraken[0]['org_name']
                for v in lKraken:
                    self.seq_dict['kraken_organisms'].append(f"{v['org_name']} ({v['tax_id']}) [{v['pct']:.1f} pct]")

            #print(self.kraken_organism)
            #SeqDict['kraken_organisms'] = merge_kraken(lKraken)

            # MLST -----------------------------
            lMLST=get_MLST_Info(self.fasta_dir,self.orgbatch_id, self.run_id)
            if len(lMLST) >0:
                self.seq_dict['mlst_scheme'] = lMLST[0]['mlst_scheme']
                self.seq_dict['mlst_seqtype'] = lMLST[0]['mlst_seqtype']
                self.seq_dict['mlst_alleles'] = lMLST[0]['mlst_alleles']

            # GTDBTK -----------------------------
            lGT=get_GTDBTK_Info(self.fasta_dir,self.orgbatch_id, self.run_id)
            if len(lGT) >0:
                self.seq_dict['gtdbtk_class'] = lGT[0]['gtdbtk_class']
                self.seq_dict['gtdbtk_fastani'] = f"{lGT[0]['gtdbtk_fastani_ref']} ({lGT[0]['gtdbtk_fastani_ani']})"

            #print(SeqDict)
            if self.seq_id:
                djIDSeq = imp_IDSeq_fromDict(self.seq_dict, self.val_log, objSeq = self.seq_id)
                if djIDSeq:
                    if djIDSeq.VALID_STATUS:
                        #print(djIDSeq.VALID_STATUS)
                        if upload:
                            djIDSeq.save(user=appuser)
                    else:
                        self.val_log.show(logTypes= ['Error'])

    #-----------------------------------------------------------------------------------
    def upload_AMR(self, Methods= ['AMR Finder'], upload=False,uploaduser=None,verbose=False):
    #-----------------------------------------------------------------------------------

        appuser = None
        if uploaduser:
            appuser = ApplicationUser.get(uploaduser)

        
        if os.path.exists(self.fasta_dir):

            # SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')
            # # Sequences -----------------------------
            # SeqDict['seq_id'] = upload_GenomeSequence(OrgBatchID, RunID, SeqDict,
            #                                 vLog,upload=upload,uploaduser=uploaduser)

            # if verbose:
            #     print(f"[WGS-AMR] {FastaDir} {OrgBatchID} {RunID} ")
            #             # Sequences -----------------------------

            # AMR Finder -----------------------------
            if 'AMR Finder' in Methods:
                lAmrFinder=get_AMRFinder_Info(self.fasta_dir,self.orgbatch_id, self.run_id)
                for row in lAmrFinder:
                    row['gene_id'] = upload_Gene(row,self.val_log,upload=upload,uploaduser=uploaduser)
                    row['seq_id'] = self.seq_id

                    djAMRgt = imp_AMRGenotype_fromDict(row,self.val_log)
                    if djAMRgt.VALID_STATUS:
                        if upload:
                            djAMRgt.save(user=appuser)
                    else:
                        self.val_log.show(logTypes= ['Error'])

            # Abricate CARD -----------------------------
            if 'Abricate card' in Methods:
                lAbCard=get_Abricate_Info(self.fasta_dir,self.orgbatch_id, self.run_id,DB='card')
                for row in lAbCard:

                    row['gene_id'] = upload_Gene(row,self.val_log,upload=upload,uploaduser=uploaduser)
                    row['seq_id'] = self.seq_id

                    djAMRgt = imp_AMRGenotype_fromDict(row,self.val_log)
                    if djAMRgt.VALID_STATUS:
                        if upload:
                            djAMRgt.save(user=appuser)
                    else:
                        self.val_log.show(logTypes= ['Error'])

            # CARD RGI -----------------------------
            if 'RGI' in Methods:
                lRGI = get_RGI_Info(self.fasta_dir,self.orgbatch_id, self.run_id)
                print(lRGI)
                for row in lRGI:
                    row['gene_id'] = upload_Gene(row,self.val_log,upload=upload,uploaduser=uploaduser)
                    row['seq_id'] = self.seq_id
                    djAMRgt = imp_AMRGenotype_fromDict(row,self.val_log)
                    if djAMRgt.VALID_STATUS:
                        if upload:
                            djAMRgt.save(user=appuser)
                    else:
                        self.val_log.show(logTypes= ['Error'])


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

# #-----------------------------------------------------------------------------------
# def gen_SeqDict(OrgBatchID,RunID,SeqType,SeqMethod,Source):
# #-----------------------------------------------------------------------------------
#     SeqDict = {
#         'seq_name'   : f"{OrgBatchID}_{RunID}",
#         'orgbatch_id': OrgBatchID,
#         'run_id'     : RunID,
#         'seq_type'   : SeqType,
#         'seq_method' : SeqMethod,
#         'source'     : Source,
#         'source_code': f"{OrgBatchID}_{RunID}",
#         'source_link': f"RDM {SeqType}: {OrgBatchID}_{RunID}",
#         'seq_file'   : 'Contigs',
#         'reference'  : ''
#     }
#     return(SeqDict)

#-----------------------------------------------------------------------------------
def upload_GenomeSequence(WGS,valLog=None, upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------
    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    djSeq = imp_Sequence_fromDict(WGS.seq_dict,vLog) 
    if djSeq.VALID_STATUS:
        if upload:
            djSeq.save(user=appuser)
        WGS.seq_id = djSeq    
    else:
        valLog.show(logTypes= ['Error'])
 
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
def upload_Trim(WGS, valLog=None, upload=False, uploaduser=None, verbose=False):  
#-----------------------------------------------------------------------------
    lstFastQC = []
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(TrimDir):

        #SeqDict = gen_SeqDict(OrgBatchID, RunID,'WGS','Illumina','CO-ADD')

        # Sequences -----------------------------
        SeqDict['seq_id'] = upload_GenomeSequence(WGS)
                                        # valLog=valLog,upload=upload,uploaduser=uploaduser)

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
def upload_CheckM(WGS, **kwargs):
#    vLog, upload=False,uploaduser=None,verbose=False):  
#-----------------------------------------------------------------------------
    upload = kwargs.get('upload',False)
    overwrite = kwargs.get('overwrite',False)
    uploaduser = kwargs.get('uploaduser',None)
    valLog = kwargs.get('valLog',None)

    lstCheckM = []

    #vLog = validation_log.Validation_Log('WGS-Assembly')

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    if os.path.exists(WGS.assembly_dir):

        # SeqDict = gen_SeqDict(WGS.orgbatch_id, WGS.run_id,'WGS',WGS.wgs_method,'CO-ADD')

        # Sequences -----------------------------
        upload_GenomeSequence(WGS,)

        if verbose:
            print(f"[WGS-Assembly] {WGS.assembly_dir} {WGS.orgbatch_id} {WGS.run_id} ")

        #sDict = {'seq_name':f"{WGS.orgbatch_id}_{WGS.run_id}"}

        # CheckM -----------------------------
        lCheckM=get_CheckM_Info(WGS)
        for row in lCheckM:
            djCheckM = imp_CheckM_fromDict(row, vLog, WGS.seq_id)
            #print(djCheckM.VALID_STATUS)
            if djCheckM.VALID_STATUS:
                
                if upload:
                    djCheckM.save(user=appuser)
                else:
                    vLog.show(logTypes= ['Error'])
            lstCheckM.append(dict(sDict,**row))
    return(lstCheckM)    


#-----------------------------------------------------------------------------------
def split_kraken(lstKraken):
#-----------------------------------------------------------------------------------
    idLst = []
    for v in lstKraken:
        _k = {}
        m = re.search(r'(.*?) \((\d+)\) \[(.*?) pct\]',v)
        if m:
            _k['org_name'] = m.group(1)
            _k['tax_id'] = int(m.group(2))
            _k['pct'] = float(m.group(3))
            idLst.append(_k)
    return(idLst)

#-----------------------------------------------------------------------------------
def agg_kraken(lstKraken,cutoff=10):
#-----------------------------------------------------------------------------------
    _agg = []
    for v in lstKraken:
        if v['pct'] >= cutoff :
            _agg.append(f"{v['org_name']} [{v['pct']:.1f} pct]")
    if len(_agg)<len(lstKraken):
        _agg.append(f"+{len(lstKraken)-len(_agg)} [<{cutoff} pct]")
    return "; ".join(_agg)
#-----------------------------------------------------------------------------------
def merge_kraken(lstKraken):
#-----------------------------------------------------------------------------------
    idLst = []
    for v in lstKraken:
        # Formatted String for ArrayField
        idLst.append(f"{v['org_name']} ({v['tax_id']}) [{v['pct']:.1f} pct]")
    return(idLst)

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
        SeqDict['kraken_organisms'] = merge_kraken(lKraken)

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

        #print(SeqDict)
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