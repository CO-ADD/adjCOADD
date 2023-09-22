import os
from datetime import datetime
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from dscreen.models import Screen_Run
from dgene.models import Genome_Sequence,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM, Gene, AMR_Genotype

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import *

# ----------------------------------------------------------------------------------------------------
def imp_Sequence_fromDict(iDict,valLog):
    """
    Create Sequence instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    validStatus = True
    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error','Organism Batch does not Exists',iDict['orgbatch_id'],'Use existing OrganismBatch ID')
        validStatus = False

    SeqType = Dictionary.get(Genome_Sequence.Choice_Dictionary["seq_type"],iDict['seq_type'])
    if SeqType is None:
        valLog.add_log('Error','Seq Type not Correct',iDict['seq_type'])
        validStatus = False

    SeqMethod = Dictionary.get(Genome_Sequence.Choice_Dictionary["seq_method"],iDict['seq_method'])
    if SeqMethod is None:
        valLog.add_log('Error','Seq Method not Correct',iDict['seq_method'])
        validStatus = False

    RunID = Screen_Run.get(iDict['run_id'])
    if RunID is None:
        valLog.add_log('Error','RunID does not Exists',iDict['run_id'])
        validStatus = False

    # Find Instance if exist
    djSeq = Genome_Sequence.get(None,iDict['seq_name'])
    if djSeq is None:
        djSeq = Genome_Sequence()
        djSeq.orgbatch_id = OrgBatch
        djSeq.seq_type = SeqType
        djSeq.seq_method = SeqMethod
        djSeq.run_id = RunID
        djSeq.seq_name = iDict['seq_name']
        valLog.add_log('Info','New Sequence',f"{iDict['seq_name']} ")

    if 'seq_date' in iDict:
        djSeq.seq_date = iDict['seq_date']
    djSeq.source = iDict['source']
    djSeq.source_code = iDict['source_code']
    djSeq.source_link = iDict['source_link']
    djSeq.reference = iDict['reference']

    djSeq.clean_Fields()
    validDict = djSeq.validate()
    if validDict:
        #validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k)

    djSeq.VALID_STATUS = validStatus

    return(djSeq)

# ----------------------------------------------------------------------------------------------------
def imp_FastQC_fromDict(iDict, valLog, objSeq = None):
    """
    Create FastQC instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    validStatus = True
    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error','Organism Batch does not Exists',iDict['orgbatch_id'],'Use existing OrganismBatch ID')
        validStatus = False

    if objSeq is None:
        objSeq = Genome_Sequence.get(None,iDict['seq_name'])
    if objSeq is None:
        valLog.add_log('Error','Sequence does not Exists',iDict['seq_name'],'Use existing Sequence')
        validStatus = False
    else:
        iDict['seq_id'] = str(objSeq)

    # Find Instance if exist
    djInst = WGS_FastQC.get(OrgBatch,iDict['seq_id'],iDict['seq'])
    if djInst is None:
        djInst = WGS_FastQC()
        djInst.orgbatch_id = OrgBatch
        djInst.seq_id = objSeq
        djInst.seq = iDict['seq']
        valLog.add_log('Info','New FastQC',f"{iDict['orgbatch_id']} {iDict['seq_id']} {iDict['seq']}",'-')

    djInst.base_stat = iDict['basic statistics']
    djInst.base_sequal = iDict['per base sequence quality']
    djInst.tile_sequal =  iDict['per tile sequence quality']
    djInst.seq_qualsc = iDict['per sequence quality scores'] 
    djInst.base_gc =  iDict['per base sequence content']
    djInst.seq_gc = iDict['per sequence gc content']
    djInst.base_N = iDict['per base n content']
    djInst.seq_len = iDict['sequence length distribution']
    djInst.seq_dup = iDict['sequence duplication levels']
    djInst.overrep = iDict['overrepresented sequences']
    djInst.adap_cont = iDict['adapter content']

    djInst.clean_Fields()
    validDict = djInst.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k,'-')

    djInst.VALID_STATUS = validStatus

    return(djInst)


# ----------------------------------------------------------------------------------------------------
def imp_CheckM_fromDict(iDict,valLog, objSeq = None):
    """
    Create CheckM instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    validStatus = True
    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error','Organism Batch does not Exists',iDict['orgbatch_id'],'Use existing OrganismBatch ID')
        validStatus = False

    if objSeq is None:
        objSeq = Genome_Sequence.get(None,iDict['seq_name'])
    if objSeq is None:
        valLog.add_log('Error','Sequence does not Exists',iDict['seq_name'],'Use existing Sequence')
        validStatus = False
    else:
        iDict['seq_id'] = str(objSeq)

    # Find Instance if exist
    djInst = WGS_CheckM.get(OrgBatch,iDict['seq_id'],iDict['assembly'])
    if djInst is None:
        djInst = WGS_CheckM()
        djInst.orgbatch_id = OrgBatch
        djInst.seq_id = objSeq
        djInst.assembly = iDict['assembly']
        valLog.add_log('Info','New CheckM',f"{iDict['orgbatch_id']} {iDict['seq_id']} {iDict['assembly']}",'-')

    djInst.assembly_qc = iDict['assembly_qc']
    djInst.marker_lineage = iDict['marker_lineage']
    djInst.n_genomes = int(iDict['n_genomes'])
    djInst.n_predit_genes =  int(iDict['n_predit_genes'])
    djInst.n_markers = int(iDict['n_markers'])
    djInst.n_marker_sets =  int(iDict['n_marker_sets'])
    djInst.genome_size = int(iDict['genome_size'])
    djInst.coding_density = round(float(iDict['coding_density']),3)
    djInst.completeness = round(float(iDict['completeness']),2)
    djInst.contamination = round(float(iDict['contamination']),2)
    djInst.gc = round(float(iDict['gc']),3)
    djInst.gc_std = round(float(iDict['gc_std']),3)
    djInst.n_ambig_bases = int(iDict['n_ambig_bases'])
    djInst.trans_table = int(iDict['trans_table'])
    djInst.n_contigs = int(iDict['n_contigs'])
    djInst.longest_contig = int(iDict['longest_contig'])
    djInst.mean_contigs = round(float(iDict['mean_contigs']),2)
    djInst.n50_contigs = int(iDict['n50_contigs'])

    if validStatus:
        djInst.clean_Fields()
        validStatus = True
        validDict = djInst.validate()
        if validDict:
            #validStatus = False
            for k in validDict:
                valLog.add_log('Warning',validDict[k],k,'-')

    djInst.VALID_STATUS = validStatus

    return(djInst)

# ----------------------------------------------------------------------------------------------------
def imp_IDSeq_fromDict(iDict,valLog, objSeq = None):
    """
    Create CheckM instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    validStatus = True
    OrgBatch = Organism_Batch.get(iDict['orgbatch_id']) 
    if OrgBatch is None:
        valLog.add_log('Error','Organism Batch does not Exists',iDict['orgbatch_id'],'Use existing OrganismBatch ID')
        validStatus = False

    if objSeq is None:
        objSeq = Genome_Sequence.get(None,iDict['seq_name'])

    if objSeq is None:
        valLog.add_log('Error','Sequence does not Exists',iDict['seq_name'],'Use existing Sequence')
        validStatus = False
    else:
        iDict['seq_id'] = str(objSeq)

    SeqFile = Dictionary.get(ID_Sequence.Choice_Dictionary["seq_file"],iDict['seq_file'])
    if SeqFile is None:
        valLog.add_log('Error','ID Type not Correct',iDict['seq_file'],'-')
        validStatus = False

    djInst = ID_Sequence.get(OrgBatch,iDict['seq_file'],iDict['seq_id'])
    if djInst is None:
        djInst = ID_Sequence()
        djInst.orgbatch_id = OrgBatch
        djInst.seq_id = objSeq
        djInst.seq_file = SeqFile
        valLog.add_log('Info','New ID-Seq',f"{iDict['orgbatch_id']} {iDict['seq_file']} ")

    djInst.kraken_organisms =iDict['kraken_organisms']
    if 'mlst_scheme' in iDict:
        djInst.mlst_scheme =iDict['mlst_scheme']
        djInst.mlst_seqtype =iDict['mlst_seqtype']
        djInst.mlst_alleles =iDict['mlst_alleles']

    if 'gtdbtk_class' in iDict:
        djInst.gtdbtk_class =iDict['gtdbtk_class']
        djInst.gtdbtk_fastani =iDict['gtdbtk_fastani']

    if 'id_notes' in iDict:
        djInst.id_notes = iDict['id_notes']
    if 'source' in iDict:
        djInst.source = iDict['source']

    if validStatus:
        djInst.clean_Fields()
        validStatus = True
        validDict = djInst.validate()
        if validDict:
            #validStatus = False
            for k in validDict:
                valLog.add_log('Warning',validDict[k],k,'-')
    djInst.VALID_STATUS = validStatus

    return(djInst)

# ----------------------------------------------------------------------------------------------------
def imp_Gene_fromDict(iDict,valLog, objSeq = None):
    """
    Create CheckM instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------

    validStatus = True
  
    if objSeq is None:
        objSeq = Genome_Sequence.get(None,iDict['seq_name'])

    if objSeq is None:
        valLog.add_log('Error','Sequence does not Exists',iDict['seq_name'],'Use existing Sequence')
        validStatus = False
    else:
        iDict['seq_id'] = str(objSeq)

    GeneType = Dictionary.get(Gene.Choice_Dictionary["gene_type"],iDict['gene_type'])
    if GeneType is None:
        valLog.add_log('Error','Gene Type not Correct',iDict['gene_type'],'-')
        validStatus = False

    # Find Instance if exist
    djGene = Gene.get(None,iDict['gene_code'])
    if djGene is None:
        djGene = Gene()
        djGene.gene_code  = iDict['gene_code']

    djGene.gene_name = iDict['gene_name']
    djGene.source = iDict['source']
 
    djGene.gene_type = GeneType
    djGene.gene_subtype = iDict['gene_subtype']
    djGene.amr_class = iDict['amr_class']
    djGene.amr_subclass = iDict['amr_subclass']

    djGene.clean_Fields()
    validDict = djGene.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k)

    djGene.VALID_STATUS = validStatus

    return(djGene)

# ----------------------------------------------------------------------------------------------------
def imp_AMRGenotype_fromDict(iDict,valLog):
    """
    Create CheckM instance from zAssembly Parser
    """
# ----------------------------------------------------------------------------------------------------
    # Find Instance if exist

    validStatus = True
    djAMRGt = AMR_Genotype.get(iDict['gene_id'],iDict['amr_method'],iDict['seq_id'],None)
    if djAMRGt is None:
        djAMRGt = AMR_Genotype()
        djAMRGt.gene_id  = iDict['gene_id']
        djAMRGt.amr_method  = iDict['amr_method']
        djAMRGt.seq_id  = iDict['seq_id']

    djAMRGt.orgbatch_id  = iDict['seq_id'].orgbatch_id
    djAMRGt.seq_coverage = round(float(iDict['amrfinder_coverage']),2)
    djAMRGt.seq_identity = round(float(iDict['amrfinder_identitiy']),2)
    djAMRGt.closest_id = iDict['amrfinder_closest']
    djAMRGt.closest_name = iDict['amrfinder_closestname']
    #djAMRGt.contig = iDict['amrfinder_contigid']

    djAMRGt.clean_Fields()
    validDict = djAMRGt.validate()
    if validDict:
        #validStatus = False
        for k in validDict:
            valLog.add_log('Warning',validDict[k],k)

    djAMRGt.VALID_STATUS = validStatus

    return(djAMRGt)
