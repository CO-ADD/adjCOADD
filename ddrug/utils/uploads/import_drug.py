import os
from pathlib import Path
from django_rdkit.models import *
#from django_rdkit.config import config
#from django.conf import settings

#from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug
#from dscreen.models import Screen_Run
from ddrug.utils.molecules import *

from apputil.models import ApplicationUser, Dictionary
from apputil.utils.data import split_StrList

# ----------------------------------------------------------------------------------------------------
def imp_Drug_fromDict(iDict,valLog):
    """
    Create Drug instance from {Dict}
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    # Remove nan
    for c in iDict:
        if iDict[c] != iDict[c]:
            iDict[c] = None

    # Find Instance if exist
    djDrug = Drug.get(iDict['drug_name'],None)
    if djDrug is None:
        djDrug = Drug()
        #djDrug.drug_id = iDict['drug_id']
        djDrug.drug_name = iDict['drug_name']
        valLog.add_log('Info',"",f"{iDict['drug_name']}",'New Drug','-') 
    djDrug.drug_type = Dictionary.get(djDrug.Choice_Dictionary["drug_type"],iDict['drug_type'])

    djDrug.n_compounds = iDict['ncmpd']
    djDrug.drug_othernames = split_StrList(iDict['drug_othernames'])
    djDrug.drug_codes = split_StrList(iDict['drug_code'])
    djDrug.drug_note = iDict['drug_note']
    djDrug.drug_panel = split_StrList(iDict['panel'])

    djDrug.antimicro = iDict['antimicro']
    djDrug.antimicro_class = iDict['antimicro_class']
    djDrug.drug_class = iDict['drug_class']
    djDrug.drug_subclass = iDict['drug_subclass']
    djDrug.drug_target = iDict['drug_target']
    #djDrug.drug_subtarget = iDict['drug_subtarget']
    #djDrug.moa = iDict['moa']

    djDrug.vendor = iDict['vendor']
    djDrug.vendor_catno = iDict['vendor_catno']

    djDrug.chembl = iDict['chembl']
    djDrug.drugbank = iDict['drugbank']
    djDrug.cas = iDict['cas']
    djDrug.pubchem = iDict['pubchem']
    djDrug.chemspider = iDict['chemspider']
    djDrug.unii = iDict['unii']
    djDrug.kegg = iDict['kegg']
    djDrug.comptox = iDict['comptox']
    djDrug.echa = iDict['echa']
    djDrug.chebi = iDict['chebi']
    djDrug.uq_imb = iDict['imb']

    djDrug.smiles = iDict['smiles']
    djDrug.smol = iDict['smiles']
    djDrug.smol = smiles2mol(iDict['smiles'],verbose=1)

    djDrug.clean_Fields()
    validStatus = True
    validDict = djDrug.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning','',k,validDict[k],'-')

    djDrug.VALID_STATUS = validStatus

    return(djDrug)

