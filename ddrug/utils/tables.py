import os
import pandas as pd
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from django.conf import settings

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST, MIC_COADD, MIC_Pub
from ddrug.utils.molecules import *

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import *

# -----------------------------------------------------------------------------------------
def get_Antibiogram(OrgID):
    """
    get MIC values for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgMIC = []

    vMIC = VITEK_AST.objects.filter(orgbatch_id__organism_id=OrgID)
    for m in vMIC:
        aDict = {}
        aDict['Drug Name'] = m.drug_name
        aDict['Drug Code'] = m.drug_codes[0]
        aDict['OrganismID'] = '_'.join(m.orgbatch_id.split('_')[0:2])
        aDict['BatchID'] = m.orgbatch_id.split('_')[2]
        if '>=' in m['MIC']:
            aDict['MIC'] = m.mic.replace('>= ','>')
        elif '<=' in m['MIC']:
            aDict['MIC'] = m.mic.replace('<= ','<=')
        else:
            aDict['MIC'] = m.mic
        aDict['BP Profile'] =m.bp_profile
        aDict['BP Source'] = m.bp_source
        aDict['Source'] = "Vitek"
        aDict['Sel Organism'] = m.organism
        aDict['Card'] = m.card_code
        if m['MIC']:
            orgMIC.append(aDict)

    pMIC = MIC_Pub.objects.filter(organism_id=OrgID)
    for m in pMIC:
        aDict = {}
        aDict['Drug Name'] = m['DRUG_NAME']
        aDict['Drug Code'] = m['DRUG_CODES']
        aDict['OrganismID'] = '_'.join(m['ORGBATCH_ID'].split('_')[0:2])
        aDict['BatchID'] = m['ORGBATCH_ID'].split('_')[2]
        aDict['MIC'] = m['MIC']
        aDict['BP Profile'] = '-'
        aDict['BP Source'] = m['RUN_ID']
        aDict['Source'] = "CO-ADD"
        aDict['Sel Organism'] = '-'
        aDict['Card'] = '-'
        if m['MIC']:
            orgMIC.append(aDict)

    cMIC = MIC_COADD.objects.filter(orgbatch_id__organism_id=OrgID)
    for m in cMIC:
        aDict = {}
        aDict['Drug Name'] = m['DRUG_NAME']
        aDict['Drug Code'] = m['DRUG_CODES']
        aDict['OrganismID'] = m['ORAGNISM_ID']
        aDict['BatchID'] = 'pub'
        if '<=' in m['MIC']:
            aDict['MIC'] = m['MIC']
        elif '<' in m['MIC']:
            aDict['MIC'] = m['MIC'].replace('<','<=')
        else:
            aDict['MIC'] = m['MIC']
        aDict['BP Profile'] =m['BP_PROFILE']
        aDict['BP Source'] = '-'
        aDict['Source'] = m['SOURCE']
        aDict['Sel Organism'] = '-'
        aDict['Card'] = '-'
        if m['MIC']:
            orgMIC.append(aDict)

    df = pd.DataFrame(orgMIC)
    pivOrg = df.pivot_table(index=["OrganismID","BatchID"],columns=["Drug Name","Source"],values=["MIC"], aggfunc=lambda x: "; ".join(list(x)))
    return(pivOrg)