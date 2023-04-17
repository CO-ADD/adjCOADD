#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

import django
from oraCastDB import oraCastDB
from zUtils import zData

from apputil.models import ApplicationUser, Dictionary
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST
from apputil.utils import slugify

from rdkit import Chem
from rdkit.Chem import AllChem
#from django_rdkit.models import *


#-----------------------------------------------------------------------------------
def reformat_OrganismID(OrgID):
#-----------------------------------------------------------------------------------
    xStr = OrgID.split("_")
    return(f"{xStr[0]}_{int(xStr[1]):04d}")

#-----------------------------------------------------------------------------------
def split_StrList(strList,sep=";"):
#-----------------------------------------------------------------------------------
    if strList:
        retLst = str(strList).split(sep)
        for i in range(len(retLst)):
            retLst[i] = retLst[i].strip()
    else:
        retLst = None
    return(retLst)

#-----------------------------------------------------------------------------------
def smiles2mol(Smiles,verbose=0):
#-----------------------------------------------------------------------------------
    try:
        xmol = Chem.MolFromSmiles(Smiles)
    except:
        xmol = None
        if verbose:
            print(f"[Invalid SMILES] {Smiles} ")
    return(xmol)

#-----------------------------------------------------------------------------------
def update_VitekAST_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------

    firstStockDate = datetime.date(2010, 1, 1)
    micSQL = """
    Select Organism_ID, Batch_ID, Card_BarCode, Card_Code,  
        Drug_ID, Drug_Name,MIC, BP_Profile, BP_Source, BP_Comment,
        Selected_Organism, Organism_Origin,
        Filename, PageNo
    From VitekAST
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[VitekID] ... ")
    micLst = OrgDB.get_dict_list(micSQL)
    nTotal = len(micLst)
    logger.info(f"[VitekID] {nTotal} ")
    OrgDB.close()

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.exists(uploaduser)

    for mic in micLst:
        orgID = reformat_OrganismID(mic['ORGANISM_ID'])
        orgClass = orgID.split("_")[0]
        batchNo = int(mic['BATCH_ID'])
        batchID = Organism_Batch.str_OrgBatchID(orgID,batchNo)

        orgBatch = Organism_Batch.exists(batchID)
        vitekCard = VITEK_Card.exists(mic['CARD_BARCODE'])
        drug = Drug.exists(mic['DRUG_NAME'],None)

        if orgBatch and vitekCard and drug:
            xx = 0
            #djMIC = VITEK_AST.exists(mic['CARD_BARCODE'])
        else:
            print(f"[NOT FOUND] {vitekCard} ({mic['CARD_BARCODE']}) {orgBatch} ({batchID}) {drug} ({mic['DRUG_NAME']}) ")



#-----------------------------------------------------------------------------------
def update_VitekID_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------

    firstStockDate = datetime.date(2010, 1, 1)
    idSQL = """
    Select Organism_ID, Batch_ID, Card_BarCode, Card_Code,  
        ID_Organism, ID_Probability, ID_Confidence,
        Filename, PageNo
    From VitekID
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[VitekID] ... ")
    IDLst = OrgDB.get_dict_list(idSQL)
    nTotal = len(IDLst)
    logger.info(f"[VitekID] {nTotal} ")
    OrgDB.close()

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.exists(uploaduser)

    for id in IDLst:
        #if card['ORGANISM_ID']:
        orgID = reformat_OrganismID(id['ORGANISM_ID'])
        orgClass = orgID.split("_")[0]
        batchNo = int(id['BATCH_ID'])
        batchID = Organism_Batch.str_OrgBatchID(orgID,batchNo)
        orgBatch = Organism_Batch.exists(batchID)
        vitekCard = VITEK_Card.exists(id['CARD_BARCODE'])

        if orgBatch and vitekCard:
            djID = VITEK_ID.exists(id['CARD_BARCODE'])
            if djID is None:
                djID = VITEK_ID()
                djID.orgbatch_id = orgBatch
                djID.card_barcode = vitekCard
            djID.id_organism = id['ID_ORGANISM']
            djID.id_probability = id['ID_PROBABILITY']
            djID.id_confidence = id['ID_CONFIDENCE']
            djID.filename = id['FILENAME']
            djID.page_no = id['PAGENO']

            djID.clean_Fields()
            validDict = djID.validate()
            if validDict:
                logger.info(f" XX {djID} {validDict} ")
            #--- Upload ---------------------------------------------------------
            if upload:
                #print(f" -> {djID} as {appuser}")
                djID.save(user=appuser)
            # else:
            #     logger.info(f" >r {djID} as {appuser}")
        else:
            print(f"[NOT FOUND] {vitekCard} ({id['CARD_BARCODE']}) {orgBatch} ({batchID}) ")

#-----------------------------------------------------------------------------------
def update_VitekCard_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------

    firstStockDate = datetime.date(2010, 1, 1)
    cardSQL = """
    Select Organism_ID, Batch_ID,  
        Card_Type, Card_Code, Card_Barcode, Expiry_Date,
        Instrument, Processing_Date, Analysis_Time
    From VitekCard
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[VitekCard] ... ")
    cardLst = OrgDB.get_dict_list(cardSQL)
    nTotal = len(cardLst)
    logger.info(f"[VitekCard] {nTotal} ")
    OrgDB.close()

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.exists(uploaduser)

    for card in cardLst:
        #if card['ORGANISM_ID']:
        orgID = reformat_OrganismID(card['ORGANISM_ID'])
        orgClass = orgID.split("_")[0]
        batchNo = int(card['BATCH_ID'])
        batchID = Organism_Batch.str_OrgBatchID(orgID,batchNo)
        orgBatch = Organism_Batch.exists(batchID)

        if orgBatch:
            djCard = VITEK_Card.exists(card['CARD_BARCODE'])
            if djCard is None:
                djCard = VITEK_Card()
                djCard.orgbatch_id = orgBatch
                djCard.card_barcode = card['CARD_BARCODE'] 
                djCard.card_code = card['CARD_CODE'] 
            djCard.card_type = Dictionary.exists(djCard.Choice_Dictionary["card_type"],card['CARD_TYPE'])
            djCard.instrument = card['INSTRUMENT']
            djCard.analysis_time = card['ANALYSIS_TIME']
            djCard.expiry_date = card['EXPIRY_DATE']
            djCard.proc_date = card['PROCESSING_DATE']

            djCard.clean_Fields()
            validDict = djCard.validate()
            if validDict:
                logger.info(f" XX {djCard} {validDict} ")
            #--- Upload ---------------------------------------------------------
            if upload:
                #print(f" -> {djTax} as {appuser}")
                djCard.save(user=appuser)
            # else:
            #     logger.info(f" >r {djCard} as {appuser}")
        else:
            print(f"[NOT FOUND] {card['CARD_BARCODE']} {orgBatch} {batchID} ")

#-----------------------------------------------------------------------------------
def update_Drug_xls(XlsFile, XlsSheet=0, upload=False, uploaduser=None, lower=True):
#-----------------------------------------------------------------------------------
    rmColumns = ['chk']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)
        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            #logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        lstRows = [{k:v for k,v in m.items()} for m in dfSheet.to_dict("rows")]
        # df -> lstDict and remove null items 
        #lstRows = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='rows')]

        # check user
        appuser = ApplicationUser.get(uploaduser)

        for row in lstRows:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in row:
                    del row[rmCol]
            
            # remove nan
            for c in row:
                if row[c] != row[c]:
                    row[c] = None

            # set instance
            djDrug = Drug.get(None,row['drug_id'])
            if djDrug is None:
                djDrug = Drug()
                djDrug.drug_id = row['drug_id']
                djDrug.drug_name = row['drug_name']
            djDrug.drug_type = Dictionary.get(djDrug.Choice_Dictionary["drug_type"],row['drug_type'])

            # Mod	Drug_OtherNames	Drug_Code		Panel		SMILES_Solv	SMILES_Salt

            djDrug.n_compounds = row['ncmpd']
            djDrug.drug_othernames = split_StrList(row['drug_othernames'])
            djDrug.drug_code = split_StrList(row['drug_code'])
            djDrug.drug_note = row['drug_note']
            djDrug.drug_panel = split_StrList(row['panel'])

            djDrug.antimicro = row['antimicro']
            djDrug.antimicro_class = row['antimicro_class']
            djDrug.drug_class = row['drug_class']
            djDrug.drug_subclass = row['drug_subclass']
            djDrug.drug_target = row['drug_target']
            #djDrug.drug_subtarget = row['drug_subtarget']
            #djDrug.moa = row['moa']

            djDrug.vendor = row['vendor']
            djDrug.vendor_catno = row['vendor_catno']

            djDrug.chembl = row['chembl']
            djDrug.drugbank = row['drugbank']
            djDrug.cas = row['cas']
            djDrug.pubchem = row['pubchem']
            djDrug.chemspider = row['chemspider']
            djDrug.unii = row['unii']
            djDrug.kegg = row['kegg']
            djDrug.comptox = row['comptox']
            djDrug.echa = row['echa']
            djDrug.chebi = row['chebi']
            djDrug.uq_imb = row['imb']

            djDrug.smiles = row['smiles']
            djDrug.smol = row['smiles']
            djDrug.smol = smiles2mol(row['smiles'],verbose=1)
            #torsionbv=TORSIONBV_FP('molecule'),
            #mfp2=MORGANBV_FP('molecule'),
            #djDrug.ffp2=FEATMORGANBV_FP(row['smiles'])
            #AllChem.GetMorganFingerprintAsBitVect(djDrug.smol)

            djDrug.clean_Fields()
            validDict = djDrug.validate()
            if validDict:
                logger.info(f" XX {djDrug} {validDict} ")
            #--- Upload ---------------------------------------------------------
            if upload:
                #print(f" -> {djTax} as {appuser}")
                djDrug.save(user=appuser)
            else:
                logger.info(f" >r {djDrug} as {appuser}")
        #Drug.update_fp()
