import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

import django
#from zUtils import zData
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from apputil.utils import validation_log
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST, MIC_COADD, MIC_Pub
from dgene.models import ID_Sequence, ID_Pub, WGS_CheckM
from dgene.utils.upload_gene import split_kraken, agg_kraken

#from ddrug.utils.import_drug import *

import ddrug.utils.vitek as Vitek

def agg_ListStr(x,unique_only=True):
    if unique_only:
        x = list(set(x))
    sx = ', '.join(['{:s}'.format(i) for i in x if i is not None])
    return(sx)



#-----------------------------------------------------------------------------------
def export_OrgBatch(OutDir):
#-----------------------------------------------------------------------------------
    AuditFields = ['astatus','acreated_at','acreated_id','aupdated_at','aupdated_id','adeleted_at','adeleted_id']

    vOrgBatch = Organism_Batch.objects.filter(astatus__gte=0).values().annotate(
#        batch_id = F("card_barcode_id__orgbatch_id__batch_id"),
#        organism_id = F("organism_id"),
        orgname = F("organism_id__organism_name"),
        strain = F("organism_id__strain_ids"),
        strain_source = F("organism_id__source"),
        strain_id = F("organism_id__strain_identification"),
        )
    logger.info(f"[OrgBatch ] {len(vOrgBatch)} ID's ")

    # for v in vOrgBatch:
    #     v['piv_value'] = f"{v['id_organism']} ({v['id_probability']})"

    dfOrgBatch = pd.DataFrame(vOrgBatch).drop(columns=AuditFields)
    # pivID = dfOrgBatch.pivot_table(index=['organism_id','orgname','batch_id'],columns=['card_code'],values=['piv_value'],aggfunc=agg_ListStr)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    xlFile = os.path.join(OutDir,"OrgBatch_Data.xlsx")
    with pd.ExcelWriter(xlFile) as writer:
        dfOrgBatch.to_excel(writer, sheet_name='OrgBatch')
        # pivID.to_excel(writer, sheet_name='ID')
        # dfAST.to_excel(writer, sheet_name='VitekAST')
        # pivAST.to_excel(writer, sheet_name='AST')

#-----------------------------------------------------------------------------------
def export_Vitek(OutDir):
#-----------------------------------------------------------------------------------

    AuditFields = ['astatus','acreated_at','acreated_id','aupdated_at','aupdated_id','adeleted_at','adeleted_id']

    vIDs = VITEK_ID.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("card_barcode_id__orgbatch_id__batch_id"),
        organism_id = F("card_barcode_id__orgbatch_id__organism_id"),
        orgname = F("card_barcode__orgbatch_id__organism_id__organism_name"),
        card_type = F("card_barcode__card_type"),
        card_code = F("card_barcode__card_code"),
        )
    logger.info(f"[Vitek    ] {len(vIDs)} ID's ")

    for v in vIDs:
        v['piv_value'] = f"{v['id_organism']} ({v['id_probability']})"

    dfID = pd.DataFrame(vIDs).drop(columns=AuditFields)
    pivID = dfID.pivot_table(index=['organism_id','orgname','batch_id'],columns=['card_code'],values=['piv_value'],aggfunc=agg_ListStr)


    vASTs = VITEK_AST.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("card_barcode_id__orgbatch_id__batch_id"),
        organism_id = F("card_barcode_id__orgbatch_id__organism_id"),
        orgname = F("card_barcode__orgbatch_id__organism_id__organism_name"),
        card_type = F("card_barcode__card_type"),
        card_code = F("card_barcode__card_code"),
        drug_name = F("drug_id__drug_name"),
        drug_codes = F("drug_id__drug_codes"),
        antimicro = F("drug_id__antimicro"),
        antimicro_class = F("drug_id__antimicro_class"),
    )
    logger.info(f"[Vitek    ] {len(vASTs)} AST's ")

    for v in vASTs:
        if v['bp_profile'] or v['mic']:
            if v['bp_profile']:  
                v['piv_value'] = f"{v['mic']} [{v['mic']}]"
            else:
                v['piv_value'] = f"{v['mic']}"
        else:
            v['piv_value'] = ''

    dfAST = pd.DataFrame(vASTs).drop(columns=AuditFields)
    pivAST = dfAST.pivot_table(index=['organism_id','orgname','batch_id'],columns=['antimicro','antimicro_class','drug_name'],values=['piv_value'],aggfunc=agg_ListStr)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    xlFile = os.path.join(OutDir,"Vitek_Data.xlsx")
    with pd.ExcelWriter(xlFile) as writer:
        dfID.to_excel(writer, sheet_name='VitekID')
        pivID.to_excel(writer, sheet_name='ID')
        dfAST.to_excel(writer, sheet_name='VitekAST')
        pivAST.to_excel(writer, sheet_name='AST')


#-----------------------------------------------------------------------------------
def export_Antibiogram(OutDir):
#-----------------------------------------------------------------------------------

    AuditFields = ['astatus','acreated_at','acreated_id','aupdated_at','aupdated_id','adeleted_at','adeleted_id','id']

    vCOADDs = MIC_COADD.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("orgbatch_id__batch_id"),
        organism_id = F("orgbatch_id__organism_id"),
        orgname = F("orgbatch_id__organism_id__organism_name"),
        drug_name = F("drug_id__drug_name"),
        drug_codes = F("drug_id__drug_codes"),
        antimicro = F("drug_id__antimicro"),
        antimicro_class = F("drug_id__antimicro_class"),
        )
    logger.info(f"[Antibiogram] {len(vCOADDs)} CO-ADD MIC's ")


    vPUBs = MIC_Pub.objects.filter(astatus__gte=0).values().annotate(
        orgname = F("organism_id__organism_name"),
        drug_name = F("drug_id__drug_name"),
        drug_codes = F("drug_id__drug_codes"),
        antimicro = F("drug_id__antimicro"),
        antimicro_class = F("drug_id__antimicro_class"),
        )
    logger.info(f"[Antibiogram] {len(vPUBs)} Pub MIC's ")

    vMICs = []
    for v in vCOADDs:
        if v['bp_profile'] or v['mic']:
            if v['bp_profile'] and v['mic']:  
                v['piv_value'] = f"{v['mic']} [{v['bp_profile']}]"
            elif v['bp_profile']:
                v['piv_value'] = f"[{v['bp_profile']}]"
            elif v['mic']:
                v['piv_value'] = f"{v['mic']}"
        else:
            v['piv_value'] = ''
        vMICs.append(v)

    for v in vPUBs:
        if v['bp_profile'] or v['mic']:
            if v['bp_profile'] and v['mic']:  
                v['piv_value'] = f"{v['mic']} [{v['bp_profile']}]"
            elif v['bp_profile']:
                v['piv_value'] = f"[{v['bp_profile']}]"
            elif v['mic']:
                v['piv_value'] = f"{v['mic']}"
        else:
            v['piv_value'] = ''
        v['batch_id'] = 'pub'
        v['organism_id'] = v['organism_id_id']

        vMICs.append(v)

    dfMIC = pd.DataFrame(vMICs).drop(columns=AuditFields)
    pivMIC = dfMIC.pivot_table(index=['organism_id','orgname','batch_id'],columns=['antimicro','antimicro_class','drug_name'],values=['piv_value'],aggfunc=agg_ListStr)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    xlFile = os.path.join(OutDir,"MIC_Data.xlsx")

    with pd.ExcelWriter(xlFile) as writer:
        dfMIC.to_excel(writer, sheet_name='CO-ADD MIC')
        pivMIC.to_excel(writer, sheet_name='MIC')
    #     #dfAST.to_excel(writer, sheet_name='Pub MIC')
    #     #pivAST.to_excel(writer, sheet_name='MICpub')

#-----------------------------------------------------------------------------------
def export_wgsFastA(OutDir):
#-----------------------------------------------------------------------------------
    AuditFields = ['astatus','acreated_at','acreated_id','aupdated_at','aupdated_id','adeleted_at','adeleted_id','id']

    vSEQs = ID_Sequence.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("seq_id__orgbatch_id__batch_id"),
        organism_id = F("seq_id__orgbatch_id__organism_id"),
        orgname = F("seq_id__orgbatch_id__organism_id__organism_name"),
        seq_name = F("seq_id__seq_name"),
        )
    logger.info(f"[wgsFastA] {len(vSEQs)} WGS FastA's ID")

    #vIDs = []
    for v in vSEQs:
        _k = split_kraken(v['kraken_organisms'])
        v['kraken2'] = agg_kraken(_k,5)
        v['mlst'] = f"{v['mlst_scheme']} ({v['mlst_seqtype']})"
        #vIDs.append(v)

    rmCols = ['kraken_organisms','mlst_seqtype','mlst_scheme','mlst_alleles','gtdbtk_fastani']
    dfID = pd.DataFrame(vSEQs).drop(columns=AuditFields+rmCols)


    
    vCMs = WGS_CheckM.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("seq_id__orgbatch_id__batch_id"),
        organism_id = F("seq_id__orgbatch_id__organism_id"),
        orgname = F("seq_id__orgbatch_id__organism_id__organism_name"),
        )
    logger.info(f"[wgsFastA] {len(vCMs)} WGS FastA's CheckM")
    for v in vCMs:
        v['piv_value'] = f"{v['marker_lineage']} [{v['assembly_qc']}]"


    dfCM = pd.DataFrame(vCMs).drop(columns=AuditFields)
    pivCM = dfCM.pivot_table(index=['organism_id','orgname','batch_id'],columns=['assembly'],values=['piv_value'],aggfunc=agg_ListStr)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)
    xlFile = os.path.join(OutDir,"Identification_Data.xlsx")

    with pd.ExcelWriter(xlFile) as writer:
        dfID.to_excel(writer, sheet_name='ID')
        #pivMIC.to_excel(writer, sheet_name='MIC')
        dfCM.to_excel(writer, sheet_name='CheckM')
        pivCM.to_excel(writer, sheet_name='Marker')

#-----------------------------------------------------------------------------------
def export_wgsAMR(OutDir):
#-----------------------------------------------------------------------------------
    pass