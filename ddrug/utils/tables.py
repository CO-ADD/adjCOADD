import os
import pandas as pd
import numpy as np

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
def get_Antibiogram_byOrgID(OrgID):
    """
    get MIC values for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgMIC = []
    showCol = ['Drug Name','Drug Class','BatchID','Source','MIC','BP Profile']
    grbyCol = ['Drug Name','Drug Class','BatchID','Source']

    OrgObj = Organism.objects.get(organism_id=OrgID)

    vMIC = VITEK_AST.objects.filter(card_barcode__orgbatch_id__organism_id=OrgObj)
    for m in vMIC:
        #print(f" -*- {m}")
        aDict = {}
        aDict['Drug Name'] = m.drug_id.drug_name
        aDict['Drug Code'] = m.drug_id.drug_codes
        aDict['Drug Class'] = m.drug_id.antimicro_class
        OrgBatchID = str(m.card_barcode.orgbatch_id)
        aDict['OrganismID'] = '_'.join(OrgBatchID.split('_')[0:2])
        aDict['BatchID'] = OrgBatchID.split('_')[2]
        if '>=' in m.mic:
            aDict['MIC'] = m.mic.replace('>= ','>')
        elif '<=' in m.mic:
            aDict['MIC'] = m.mic.replace('<= ','<=')
        else:
            aDict['MIC'] = m.mic
        aDict['BP Profile'] =m.bp_profile
        aDict['BP Source'] = m.bp_source
        aDict['Source'] = "Vitek"
        aDict['Sel Organism'] = m.organism
        aDict['Card'] = m.card_barcode.card_code
        orgMIC.append(aDict)

    pMIC = MIC_Pub.objects.filter(organism_id=OrgObj)
    for m in pMIC:
        aDict = {}
        aDict['Drug Name'] = m.drug_id.drug_name
        aDict['Drug Code'] = m.drug_id.drug_codes
        aDict['Drug Class'] = m.drug_id.antimicro_class
        aDict['OrganismID'] = OrgID
        aDict['BatchID'] = 'pub'
        aDict['MIC'] = m.mic
        aDict['BP Profile'] = m.bp_profile
        aDict['BP Source'] = '-'
        aDict['Source'] = m.source
        aDict['Sel Organism'] = '-'
        aDict['Card'] = '-'
        orgMIC.append(aDict)

    cMIC = MIC_COADD.objects.filter(orgbatch_id__organism_id=OrgObj)
    for m in cMIC:
        aDict = {}
        aDict['Drug Name'] = m.drug_id.drug_name
        aDict['Drug Code'] = m.drug_id.drug_codes
        aDict['Drug Class'] = m.drug_id.antimicro_class
        OrgBatchID = str(m.orgbatch_id)
        aDict['OrganismID'] = '_'.join(OrgBatchID.split('_')[0:2])
        aDict['BatchID'] = OrgBatchID.split('_')[2]
        if '<=' in m.mic:
            aDict['MIC'] = m.mic
        elif '<' in m.mic:
            aDict['MIC'] = m.mic.replace('<','<=')
        else:
            aDict['MIC'] = m.mic
        aDict['BP Profile'] =m.bp_profile
        aDict['BP Source'] = m.run_id
        aDict['Source'] = "CO-ADD"
        aDict['Sel Organism'] = '-'
        aDict['Card'] = '-'
        orgMIC.append(aDict)

    df = pd.DataFrame(orgMIC)
    #df.to_excel(f"{OrgID}_Antibio.xlsx")
    df = df.fillna("-")

    agg_df = df[showCol].groupby(grbyCol).aggregate(lambda x: ", ".join(list(np.unique(x))))

    return(agg_df)
