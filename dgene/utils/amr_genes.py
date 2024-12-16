import os
import pandas as pd
import numpy as np

from pathlib import Path
from django.http import HttpResponse

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from dgene.models import AMR_Genotype

# -----------------------------------------------------------------------------------------
def get_AMRGenes_byOrgID_Html(pk, with_style = False):

    displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source']

    print(f"[get_AMRGenes_byOrgID] {pk}")
    df = get_AMRGenes_byOrgID(str(pk))
    if df is not None:
        df.reset_index(inplace=True)
        #df = df[displaycols]
        df_entries=len(df)

        print(f"[get_AMRGenes_byOrgID] {pk} : {df_entries} ")
        piv_table = piv_AMRGenes_byOrgID(df)
        print(f"[get_AMRGenes_byOrgID] {pk} : {df_entries} -> {len(piv_table)}")

        if with_style:
            html_table=df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
            html_pivtable = piv_table.to_html()
        else:
            html_table=df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
            html_pivtable = piv_table.to_html()

        table={'n_entries':df_entries, 'html_table': html_table, 'pivot_table': html_pivtable} 
    else:
        table = {'n_entries':0, 'html_table': None, 'pivot_table': None}
    return table 

# -----------------------------------------------------------------------------------------
def Export_AMRGenes_byOrgID_byOrgID(request, pk):
    displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source']
    xlsx_name = f"AMRGenes_{str(pk)}.xlsx"

    print(f"[get_AMRGenes_byOrgID] {pk}")
    df = get_AMRGenes_byOrgID(str(pk))
    if df is not None:
        df.reset_index(inplace=True)
        df = df[displaycols]
        df_entries=len(df)
        
        print(f"[get_AMRGenes_byOrgID] {pk} : {df_entries} ")
        piv_table = piv_AMRGenes_byOrgID(df)
        print(f"[get_AMRGenes_byOrgID] {pk} : {df_entries} -> {len(piv_table)} -> {xlsx_name}")

        response = HttpResponse(content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = f'attachment; filename={xlsx_name}'

        with pd.ExcelWriter(response) as writer:
            df.to_excel(writer, sheet_name='Data')
            piv_table.to_excel(writer, sheet_name='Pivot')

    return response

# -----------------------------------------------------------------------------------------
def piv_AMRGenes_byOrgID(df):
    piv_table = df.pivot_table(index=['Gene SubType', 'AMR Class'], columns=['BatchID','Method'], values=['Gene Name'],  
                                aggfunc= lambda x:  "; ".join([str(y) for y in x]))
    piv_table = piv_table.fillna("-").astype(str)
    return(piv_table)

# -----------------------------------------------------------------------------------------
def get_AMRGenes_byOrgID(OrgID):
    """
    get AMR Genes for Organism ID and prepare aggregate table as Dataframe 
    """
# -----------------------------------------------------------------------------------------
    orgGene = []
    OrgObj = Organism.objects.get(organism_id=OrgID)

    vAMR = AMR_Genotype.objects.filter(seq_id__orgbatch_id__organism_id=OrgObj)
    #print(f"Getting {len(vMIC)} Vitek AST data for {OrgID} ")
    for m in vAMR:
        aDict = {}
        OrgBatchID = str(m.seq_id.orgbatch_id)
        aDict['OrganismID'] = '_'.join(OrgBatchID.split('_')[0:2])
        aDict['BatchID'] = OrgBatchID.split('_')[2]
        aDict['Gene Name'] = m.gene_id.gene_code
        aDict['Gene SubType'] = m.gene_id.gene_subtype
        aDict['AMR Class'] = m.gene_id.amr_class
        aDict['Method'] = m.amr_method
        aDict['Coverage'] = m.seq_coverage
        aDict['Identity'] = m.seq_identity
        orgGene.append(aDict)


    if len(orgGene) > 0:
        df = pd.DataFrame(orgGene)
        df = df.fillna("-").astype(str)
        return(df)
    else:
        print(" No AMR Gene data found")
        return(None)
