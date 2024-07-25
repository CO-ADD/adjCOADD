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
from ddrug.utils.bio_data import agg_DR

#from apputil.utils.data import *
from apputil.utils.data_style import highlight_RSI

# -----------------------------------------------------------------------------------------
def get_Antibiogram_byOrgID_Html(pk, displaycols, with_style = False):

    df = get_Antibiogram_byOrgID(str(pk))
    if df is not None:
        df.reset_index(inplace=True)
        df = df[displaycols]
        df_entries=len(df)
        
        piv_table = piv_Antibiogram_byOrgID(df)

        # Styling pivottable
        #print(f"HMTL {len(piv_table)} ")

        if with_style:       
            html_table=df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
            
            html_pivtable = piv_table.style.applymap(highlight_RSI)
            html_pivtable = html_pivtable.set_table_attributes('class="table table-bordered fixTableHead"').to_html()
        else:
            html_table    = df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
            html_pivtable = piv_table.to_html()

        table={'n_entries':df_entries, 'html_table': html_table, 'pivot_table': html_pivtable} 
    else:
        table = {'n_entries':0, 'html_table': None, 'pivot_table': None}
    return table 

# -----------------------------------------------------------------------------------------
def piv_Antibiogram_byOrgID(df):
    print(f"Pivot {len(df)} ")

    piv_table = df.pivot_table(columns='BatchID',index=['Drug Class', 'Drug Name', ], values=['BP Profile', 'MIC'],  
                                aggfunc= lambda x:  " ".join([str(y) for y in x]))
    #.sort_values(by=['Drug Class'],ascending=False)
    piv_table = piv_table.fillna("-").astype(str)
    return(piv_table)
    

# -----------------------------------------------------------------------------------------
def get_Antibiogram_byOrgID(OrgID):
    """
    get MIC values for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgMIC = []
    # showCol = ['Drug Class','Drug Name','BatchID','Source','MIC','Profile','BP Source']
    # grbyCol = ['Drug Name','Drug Class','BatchID','Source']

    OrgObj = Organism.objects.get(organism_id=OrgID)

    vMIC = VITEK_AST.objects.filter(card_barcode__orgbatch_id__organism_id=OrgObj)
    #print(f"Getting {len(vMIC)} Vitek AST data for {OrgID} ")
    for m in vMIC:
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
        aDict['BP Source'] = m.card_barcode.card_code
        aDict['Source'] = "Vitek"
        orgMIC.append(aDict)

    pMIC = MIC_Pub.objects.filter(organism_id=OrgObj)
    #print(f"Getting {len(pMIC)} MIC Pub data for {OrgID} ")
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
        orgMIC.append(aDict)

    cMIC = MIC_COADD.objects.filter(orgbatch_id__organism_id=OrgObj)
    #print(f"Getting {len(cMIC)} MIC CO-ADD data for {OrgID} ")
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
        orgMIC.append(aDict)
        #print(aDict)
        
    if len(orgMIC) > 0:
        #print(f"GroupBy {len(orgMIC)} Dataframe for {OrgID} ")
        df = pd.DataFrame(orgMIC)
        #df.to_excel(f"{OrgID}_Antibio.xlsx")
        df = df.fillna("-").astype(str)

        showCol = ['Drug Name','Drug Class','BatchID','Source','MIC','BP Profile','BP Source']
        grbyCol = ['Drug Name','Drug Class','BatchID','Source']

        agg_df = df[showCol].groupby(grbyCol) \
                            .aggregate(lambda x: ", ".join(list(np.unique(x)))).sort_values(by=['Drug Class'],ascending=True)
                            # .aggregate(lambda x: agg_DR(x))                 
        return(agg_df)
    else:
        print(" No MIC data found")
        return(None)


# -----------------------------------------------------------------------------------------
def get_Identification_byOrgID(OrgID):
    """
    get Identification values for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgMIC = []
    showCol = ['BatchID','Source','Organism Name','Method','Date']
    grbyCol = ['BatchID','Source']

    OrgObj = Organism.objects.get(organism_id=OrgID)

    #vMIC = VITEK_AST.objects.filter(card_barcode__orgbatch_id__organism_id=OrgObj)

# -----------------------------------------------------------------------------------------
def get_Genes_byOrgID(OrgID):
    """
    get MIC values for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgMIC = []
    showCol = ['Gene Name','Gene Class','BatchID','Source','MIC','BP Profile','BP Source']
    grbyCol = ['Drug Name','Drug Class','BatchID','Source']

    OrgObj = Organism.objects.get(organism_id=OrgID)

