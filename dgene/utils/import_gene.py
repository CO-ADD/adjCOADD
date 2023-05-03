import os
from datetime import datetime
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from dgene.models import Gene,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import *

# ----------------------------------------------------------------------------------------------------
def imp_FastQC_fromDict(iDict,valLog):
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

    # Find Instance if exist
    djInst = WGS_FastQC.get(OrgBatch,iDict['seq'],iDict['run_id'])
    if djInst is None:
        djInst = WGS_FastQC()
        djInst.orgbatch_id = OrgBatch
        djInst.seq = iDict['seq']
        djInst.run_id = iDict['run_id']
        valLog.add_log('Info','New FastQC',f"{iDict['orgbatch_id']} {iDict['seq']} {iDict['run_id']}",'-')

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
def imp_CheckM_fromDict(iDict,valLog):
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

    # Find Instance if exist
    djInst = WGS_CheckM.get(OrgBatch,iDict['run_id'])
    if djInst is None:
        djInst = WGS_CheckM()
        djInst.orgbatch_id = OrgBatch
        djInst.run_id = iDict['run_id']
        valLog.add_log('Info','New CheckM',f"{djInst.orgbatch_id} {djInst.run_id}",'-')

    djInst.marker_lineage = iDict['marker_lineage']
    djInst.n_genomes = iDict['n_genomes']
    djInst.n_predit_genes =  iDict['n_predit_genes']
    djInst.n_markers = iDict['n_markers'] 
    djInst.n_marker_sets =  iDict['n_marker_sets']
    djInst.n_scaffolds = iDict['n_scaffolds']
    djInst.genome_size = iDict['genome_size']
    djInst.coding_density = iDict['coding_density']
    djInst.completeness = iDict['completeness']
    djInst.contamination = iDict['contamination']
    djInst.gc = iDict['gc']
    djInst.gc_std = iDict['gc_std']
    djInst.n_ambig_bases = iDict['n_ambig_bases']
    djInst.longest_scaffold = iDict['longest_scaffold']
    djInst.mean_scaffols = iDict['mean_scaffols']
    djInst.n50 = iDict['n50']
    djInst.trans_table = iDict['trans_table']

    if validStatus:
        djInst.clean_Fields()
        validStatus = True
        validDict = djInst.validate()
        if validDict:
            validStatus = False
            for k in validDict:
                valLog.add_log('Warning',validDict[k],k,'-')

    djInst.VALID_STATUS = validStatus

    return(djInst)

# ----------------------------------------------------------------------------------------------------
def imp_IDSeq_fromDict(iDict,valLog):
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

    IDType = Dictionary.get(ID_Sequence.Choice_Dictionary["id_type"],iDict['id_type'])
    if IDType is None:
        valLog.add_log('Error','ID Type not Correct',iDict['id_type'],'-')
        validStatus = False

    # Find Instance if exist
    djInst = ID_Sequence.get(OrgBatch,iDict['id_type'],iDict['run_id'])
    if djInst is None:
        djInst = ID_Sequence()
        djInst.orgbatch_id = OrgBatch
        djInst.id_type = IDType
        djInst.run_id = iDict['run_id']
        valLog.add_log('Info','New ID-Seq',f"{djInst.orgbatch_id} {djInst.id_type} {djInst.run_id}",'-')

    djInst.id_method = iDict['id_method']
    djInst.id_organisms =iDict['id_organisms']
    if 'id_date' in iDict:
        djInst.id_date = iDict['id_date']
    else:
        djInst.id_date = datetime.now()
    if 'id_notes' in iDict:
        djInst.id_notes = iDict['id_notes']
    if 'source' in iDict:
        djInst.source = iDict['source']

    if validStatus:
        djInst.clean_Fields()
        validStatus = True
        validDict = djInst.validate()
        if validDict:
            validStatus = False
            for k in validDict:
                valLog.add_log('Warning',validDict[k],k,'-')

    djInst.VALID_STATUS = validStatus

    return(djInst)