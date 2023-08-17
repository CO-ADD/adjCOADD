#
#
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
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock, OrgBatch_Image
from django.utils.text import slugify

#-----------------------------------------------------------------------------------
def reformat_OrganismID(OrgID):
#-----------------------------------------------------------------------------------
    xStr = OrgID.split("_")
    return(f"{xStr[0]}_{int(xStr[1]):04d}")

#-----------------------------------------------------------------------------------
def split_StrList(strList,sep=";"):
    if strList:
        retLst = strList.split(sep)
        for i in range(len(retLst)):
            retLst[i] = retLst[i].strip()
    else:
        retLst = None
    return(retLst)

#--------------------------------------------------------------------------------------
def get_subdir(OrgBatchID,binsize=100):
    """
     Gets the SubFolder name based on the XX_NNNN with splits into 200
      GN_0000, GN_0200, GN_0400, ... ,GN_1200, GN_1400, GN_1600
    """
    _org = OrgBatchID.split('_')
    return(f"{_org[0]}_{int(int(_org[1])/binsize)*binsize:04d}")

#-----------------------------------------------------------------------------------
def delete_OrgCulture(upload=False):
#-----------------------------------------------------------------------------------
    if upload:
        logger.info(f"[adjCOADD] Remove all Org_Culture entries ")
        Organism_Culture.objects.all().delete()


def update_OrgBatchImg(JpegFolder,upload=False,uploaduser=None):
    nJpeg = 0
    if JpegFolder:
        DirFiles = os.listdir(JpegFolder)
        JpegFiles = [f for f in DirFiles if f.endswith(".jpeg")]
        nJpeg = len(JpegFiles)

    if nJpeg>0:
        vImages = []

        # check user
        appuser = ApplicationUser.get(uploaduser)

        for i in range(nJpeg):
            djImg = OrgBatch_Image.get(JpegFiles[i])
            if djImg is None:
                djImg = OrgBatch_Image()
                djImg.image_name = JpegFiles[i]
                djImg.image_type = 'image/jpeg'
                djImg.image_source = 'IMB' 
                djImg.image_desc = 'Initial agarplate'
                djImg.orgbatch_id = Organism_Batch.get(JpegFiles[i].replace(".jpeg",""))
                djImg.image_file = f"images/orgbatch/{get_subdir(JpegFiles[i])}/{JpegFiles[i]}"
            else:
                djImg.image_type = 'image/jpeg'
                djImg.image_source = 'IMB' 
                djImg.image_desc = 'Initial agarplate'
                djImg.orgbatch_id = Organism_Batch.get(JpegFiles[i].replace(".jpeg",""))
                djImg.image_file = f"images/orgbatch/{get_subdir(JpegFiles[i])}/{JpegFiles[i]}"

            djImg.clean_Fields()
            validDict = djImg.validate()

            if validDict:
                logger.info(f" XX {djImg} {validDict} ")
                #nProc['notValid'] = nProc['notValid'] + 1
            else:
                # --- Upload ---------------------------------------------------------
                if upload:
                    # if nProcessed%OutputN == 0:
                    #     eTime,sTime = nTime.remains(nProcessed)
                    #     logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} ----> {repr(djCult)}") 

                    djImg.save(user=appuser)
                    #nProc['Saved'] = nProc['Saved'] + 1 

#-----------------------------------------------------------------------------------
def update_OrgCulture_xls(XlsFile, XlsSheet=0,upload=False,uploaduser=None, lower=False, OutputN=50):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProc['notFound'] = 0


    rmColumns = ['organism_id','Org ID','Batch ID','X1Date','X2Date','Passages']
    rmColumns = []
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        cultLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(cultLst)
        nTime  = zData.Timer(nTotal)
        nProcessed = 0


        # check user
        appuser = ApplicationUser.get(uploaduser)

        for entry in cultLst:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]


        for cult in cultLst:
            nProcessed = nProcessed + 1
            #print(stock)
            djCult = Organism_Culture()
            OrganismID = cult['organism_id']
            orgID = Organism.get(OrganismID)

            if orgID is not None:
                djCult.organism_id = orgID
                cult.pop('organism_id')

                if 'biologist' in cult:
                    djCult.biologist = ApplicationUser.get(cult['biologist'])
                    cult.pop('biologist')
                if 'culture_type' in cult:
                    djCult.culture_type = Dictionary.get(djCult.Choice_Dictionary["culture_type"],cult['culture_type'])
                    cult.pop('culture_type')
                if 'culture_source' in cult:
                    djCult.culture_source = Dictionary.get(djCult.Choice_Dictionary["culture_source"],cult['culture_source'])
                    cult.pop('culture_source')

                # set values in instance
                for e in cult:
                    setattr(djCult,e,cult[e])


                djCult.clean_Fields()
                validDict = djCult.validate()

                if validDict:
                    logger.info(f" XX {djCult} {validDict} ")
                    nProc['notValid'] = nProc['notValid'] + 1
                else:
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} ----> {repr(djCult)}") 

                        djCult.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    else:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >chk> {repr(djCult)}")
            else:
                logger.info(f"[OrgCulture] Organism not Found '{OrganismID}' ")
                nProc['notFound'] = nProc['notFound'] + 1
    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[OrgCulture] [{nTotal-(nProc['Saved']+nProc['notValid']+nProc['notFound'])}] {nTotal} {sTime} {nProc}")


#-----------------------------------------------------------------------------------
def delete_OrgBatchStock(upload=False):
#-----------------------------------------------------------------------------------
    if upload:
        logger.info(f"[adjCOADD] Remove all OrgBatch_Stock entries ")
        OrgBatch_Stock.objects.all().delete()
    
#-----------------------------------------------------------------------------------
def update_OrgBatchStock_xls(XlsFile, XlsSheet=0,upload=False,uploaduser=None, lower=False, OutputN=50):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProc['notFound'] = 0


    rmColumns = ['organism_id','Org ID','Batch ID','X1Date','X2Date','Passages']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        stockLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(stockLst)
        nTime  = zData.Timer(nTotal)
        nProcessed = 0


        # check user
        appuser = ApplicationUser.get(uploaduser)

        for entry in stockLst:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]


        for stock in stockLst:
            nProcessed = nProcessed + 1
            #print(stock)
            djStock = OrgBatch_Stock()
            OrgbatchID = stock['orgbatch_id']
            orgBatch = Organism_Batch.get(OrgbatchID)

            if orgBatch is not None:
                djStock.orgbatch_id = orgBatch
                stock.pop('orgbatch_id')

                if 'biologist' in stock:
                    djStock.biologist = ApplicationUser.get(stock['biologist'])
                    stock.pop('biologist')
                if 'stock_type' in stock:
                    djStock.stock_type = Dictionary.get(djStock.Choice_Dictionary["stock_type"],stock['stock_type'])
                    stock.pop('stock_type')

                # set values in instance
                for e in stock:
                    setattr(djStock,e,stock[e])


                djStock.clean_Fields()
                validDict = djStock.validate()

                if validDict:
                    logger.info(f" XX {djStock} {validDict} ")
                    nProc['notValid'] = nProc['notValid'] + 1
                else:
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} ----> {repr(djStock)}") 

                        djStock.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    else:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >chk> {repr(djStock)}")
            else:
                logger.info(f"[OrgBatchStock] Organism Batch not Found '{OrgbatchID}' ")
                nProc['notFound'] = nProc['notFound'] + 1
    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[OrgBatchStock] [{nTotal-(nProc['Saved']+nProc['notValid']+nProc['notFound'])}] {nTotal} {sTime} {nProc}")


#-----------------------------------------------------------------------------------
def update_OrgBatchStock_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------
    fix_StockType = {
        "Master"    : "MasterStock",
        "Master_LN2":"MasterLN2",
        "Stock"     :"Stock",
        "HighUse"   :"HighUse",
    }

    firstStockDate = datetime.date(2010, 1, 1)
    orgSQL = """
    Select Organism_ID, Batch_ID,  
        Stock_ID, Stock_Date, Stock_Description, Stock_Type,
        Stock_Loc_Freezer, Stock_Loc_Column, Stock_Loc_Rack, Stock_Loc_Slot,
        Stock_nCreated, Stock_nLeft, Stock_Passagen,      
        Biologist
    From Organism_Stock
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[OrgBatch] ... ")
    batchLst = OrgDB.get_dict_list(orgSQL)
    nTotal = len(batchLst)
    logger.info(f"[OrgBatch] {nTotal} ")
    OrgDB.close()

    nProc = {}
    nProc['Saved'] = 0
    nProc['Empty'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    for batch in batchLst:
        if batch['ORGANISM_ID']:
            orgID = reformat_OrganismID(batch['ORGANISM_ID'])
            orgClass = orgID.split("_")[0]
            batchID = Organism_Batch.str_BatchID(int(batch['BATCH_ID']))
            OrgbatchID = Organism_Batch.str_OrgBatchID(orgID,batchID)

            if orgClass in ['GN','GP','FG','MB']:

                #djStock = OrgBatch_Stock.get(OrgbatchID,batch['STOCK_TYPE'])
                #iforaCastDB djStock is None:
                djStock = OrgBatch_Stock()
                orgBatch = Organism_Batch.get(OrgbatchID)
                
                #print(OrgbatchID)
                if orgBatch is not None:
                    if batch['STOCK_DATE'] is None:
                        batch['STOCK_DATE'] = firstStockDate
                    if batch['STOCK_TYPE'] in fix_StockType:
                        sType = fix_StockType[batch['STOCK_TYPE']]
                        batch['STOCK_TYPE'] = Dictionary.get(djStock.Choice_Dictionary["stock_type"],sType)
                    else:
                        logger.info(f"{batch} WRONG StockType {batch['STOCK_TYPE']} ")

                    if batch['STOCK_NLEFT'] is None:
                        batch['STOCK_NLEFT'] = 0
                    if batch['STOCK_NLEFT'] < 1:
                        nProc['Empty'] = nProc['Empty'] + 1

                    djStock = OrgBatch_Stock.get(batch['STOCK_ID'])
                    if djStock is None:
                        djStock = OrgBatch_Stock()
                        djStock.orgbatch_id = orgBatch
                        djStock.stock_date = batch['STOCK_DATE'] 
                        djStock.stock_type = batch['STOCK_TYPE']
                        djStock.stock_id = batch['STOCK_ID']

                    djStock.stock_notes = batch['STOCK_DESCRIPTION']
                    djStock.location_freezer = batch['STOCK_LOC_FREEZER']
                    djStock.location_rack = batch['STOCK_LOC_RACK']
                    djStock.location_column = batch['STOCK_LOC_COLUMN']
                    djStock.location_slot = batch['STOCK_LOC_SLOT']
                    djStock.passage_no = batch['STOCK_PASSAGEN']
                    djStock.n_created = batch['STOCK_NCREATED']
                    djStock.n_left = batch['STOCK_NLEFT']
                    djStock.biologist = ApplicationUser.get(batch['BIOLOGIST'])

                    djStock.clean_Fields()
                    validDict = djStock.validate()
                    if validDict:
                        logger.info(f"{djStock} {validDict} ")
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        logger.info(f" -> {djStock} ")
                        djStock.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    else:
                        print(f" >r {djStock}")
                        #logger.info(f" >r {djStock}")
                else:
                    logger.info(f" XX OrgBatch NOT found {batch['STOCK_ID']} {batchID} ")

            else:
                logger.info(f" XX NOT in Class {batch['STOCK_ID']} {batch['ORGANISM_ID']} ")
                nProc['notClass'] = nProc['notClass'] + 1
        else:
            logger.info(f" XX Organism NOT found {batch['STOCK_ID']} {batch['ORGANISM_ID']} ")
            nProc['notFound'] = nProc['notFound'] + 1
    logger.info(f"[OrgBatchStock] [{nTotal-(nProc['Saved']+nProc['notClass']+nProc['notFound'])}] {nTotal} {nProc}")

#-----------------------------------------------------------------------------------
def update_OrgBatch_xls(XlsFile, XlsSheet=0,upload=False,uploaduser=None, lower=False, OutputN=20):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0
    nProc['notFound'] = 0


    rmColumns = ['chk']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        batchLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(batchLst)
        nTime  = zData.Timer(nTotal)
        nProcessed = 0

        # check user
        appuser = ApplicationUser.get(uploaduser)

        for entry in batchLst:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]


        for batch in batchLst:
            nProcessed = nProcessed + 1
            OrgbatchID = batch['orgbatch_id']
            BatchID = batch['batch_id']
            OrgID = batch['organism_id']

            djOrg = Organism.get(OrgID)
            if djOrg:
                djBatch = Organism_Batch.get(OrgbatchID)
                if djBatch is None:
                    djBatch = Organism_Batch()
                    djBatch.orgbatch_id = OrgbatchID
                    djBatch.batch_id = BatchID
                    djBatch.organism_id = djOrg
                batch.pop('orgbatch_id')
                batch.pop('batch_id')
                batch.pop('organism_id')

                if 'biologist' in batch:
                    djBatch.biologist = ApplicationUser.get(batch['biologist'])
                    batch.pop('biologist')
                if 'qc_status' in batch:
                    djBatch.qc_status = Dictionary.get(djBatch.Choice_Dictionary["qc_status"],batch['qc_status'],None)
                    batch.pop('qc_status')

                # set values in instance
                for e in batch:
                    setattr(djBatch,e,batch[e])


                djBatch.clean_Fields()
                validDict = djBatch.validate()
                if validDict:
                    logger.info(f" XX {djBatch} {validDict} ")
                    nProc['notValid'] = nProc['notValid'] + 1
                else:
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} ----> {repr(djBatch)}") 

                        djBatch.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    else:
                        if nProcessed%OutputN == 0:
                            eTime,sTime = nTime.remains(nProcessed)
                            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >chk> {repr(djBatch)}")
            else:
                logger.info(f"[OrgBatch] Organism not Found '{OrgID}' ")
                nProc['notFound'] = nProc['notFound'] + 1
            
    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[OrgBatch] [{nTotal-(nProc['Saved']+nProc['notValid']+nProc['notFound'])}] {nTotal} {sTime} {nProc}")


#-----------------------------------------------------------------------------------
def update_OrgBatch_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------
    orgSQL = """
    Select Organism_ID, Batch_ID, Batch_Description, Stock_Date,
        QC, QC_Record,
        Master_Level, Stock_LEvel, Highuse_Level,      
        Biologist
    From Organism_Batch
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[OrgBatch] ... ")
    batchLst = OrgDB.get_dict_list(orgSQL)
    nTotal = len(batchLst)
    logger.info(f"[OrgBatch] {nTotal} ")
    OrgDB.close()

    # check user
    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    nProc = {}
    nProc['Saved'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0

    for batch in batchLst:
        orgID = reformat_OrganismID(batch['ORGANISM_ID'])
        orgClass = orgID.split("_")[0]
        batchID = Organism_Batch.str_BatchID(int(batch['BATCH_ID']))
        OrgbatchID = Organism_Batch.str_OrgBatchID(orgID,batchID)

        if orgClass in ['GN','GP','FG','MB']:
            djOrg = Organism.get(orgID)
            if djOrg:
                djBatch = Organism_Batch.get(OrgbatchID)
                if djBatch is None:
                    djBatch = Organism_Batch()
                    djBatch.orgbatch_id = OrgbatchID
                    djBatch.batch_id = batchID
                    djBatch.organism_id = djOrg
                djBatch.batch_notes = batch['BATCH_DESCRIPTION']
                djBatch.qc_status = Dictionary.get("QC_Status",batch['QC'])
                djBatch.qc_record = batch['QC_RECORD']
                djBatch.stock_date = batch['STOCK_DATE']
                if batch['MASTER_LEVEL'] is None:
                    batch['MASTER_LEVEL'] = 0
                if batch['STOCK_LEVEL'] is None:
                    batch['STOCK_LEVEL'] = 0
                if batch['HIGHUSE_LEVEL'] is None:
                    batch['HIGHUSE_LEVEL'] = 0
                djBatch.stock_level = [batch['MASTER_LEVEL'] ,batch['STOCK_LEVEL'] ,batch['HIGHUSE_LEVEL'] ]
                djBatch.biologist = ApplicationUser.get(batch['BIOLOGIST'])

                djBatch.clean_Fields()
                validDict = djBatch.validate()
                if validDict:
                    logger.info(f"{djBatch} {validDict} ")

                # --- Upload ---------------------------------------------------------
                if upload:
                    logger.info(f" -> {djBatch} as {appuser}")
                    djBatch.save(user=appuser)
                    nProc['Saved'] = nProc['Saved'] + 1
                else:
                    logger.info(f" >r {djBatch} as {appuser}")
            else:
                logger.info(f"[OrganismID] not Found '{orgID}' ")
                nProc['notFound'] = nProc['notFound'] + 1
        else:
            nProc['notClass'] = nProc['notClass'] + 1
    logger.info(f"[OrgBatch] [{nTotal-(nProc['Saved']+nProc['notClass']+nProc['notFound'])}] {nTotal} {nProc}")

#-----------------------------------------------------------------------------------
def update_Organism_xls(XlsFile, XlsSheet=0,upload=False,uploaduser=None, lower=False, OutputN=20):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0


    rmColumns = ['chk']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        orgLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(orgLst)

        # check user
        appuser = ApplicationUser.get(uploaduser)

        for entry in orgLst:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]

        for org in orgLst:
            orgID = org['organism_id']
            orgClass = orgID.split("_")[0]

            if orgClass in ['GN','GP','FG','MB']:
            # check if instance exists in DB
                djOrg = Organism.get(orgID)
                if djOrg is None:
                    djOrg = Organism()
                    djOrg.organism_id = orgID

                djOrg.organism_name = Taxonomy.get(org['organism_name'],verbose=1)
                org.pop('organism_name')

# strain_code	strain_ids	strain_panel	strain_notes	res_property	gen_property	sero_clone	strain_identification	
# strain_origin	source	source_code	tax_id	reference	mta_document	mta_notes	
# collect_date	collect_region	collect_country	collect_site	collect_specie	collect_tissue	patient_diagnosis	patient	
# 			

                if 'lab_restriction' in org:
                    djOrg.lab_restriction = Dictionary.get(djOrg.Choice_Dictionary["lab_restriction"],org['lab_restriction'],None)
                    org.pop('lab_restriction')
                if 'risk_group' in org:
                    djOrg.risk_group = Dictionary.get(djOrg.Choice_Dictionary["risk_group"],org['risk_group'],None)
                    org.pop('risk_group')
                if 'pathogen_group' in org:
                    djOrg.pathogen_group = Dictionary.get(djOrg.Choice_Dictionary["pathogen_group"],None,org['pathogen_group'])
                    org.pop('pathogen_group')
                if 'oxygen_pref' in org:
                    djOrg.oxygen_pref = Dictionary.get(djOrg.Choice_Dictionary["oxygen_pref"],org['oxygen_pref'],None)
                    org.pop('oxygen_pref')
                if 'biologist' in org:
                    djOrg.biologist = ApplicationUser.get(org['biologist'])
                    org.pop('biologist')
                if 'mta_status' in org:
                    djOrg.mta_status = Dictionary.get(djOrg.Choice_Dictionary['mta_status'],org['mta_status'],None)
                    org.pop('mta_status')
                if 'strain_type' in org:
                    djOrg.strain_type = Dictionary.get_DictValues_fromStrList(djOrg.Choice_Dictionary["strain_type"],org['strain_type'],None)
                    org.pop('strain_type')
                if 'strain_panel' in org:
                    djOrg.strain_panel = Dictionary.get_DictValues_fromStrList(djOrg.Choice_Dictionary["strain_panel"],org['strain_panel'],None)
                    org.pop('strain_panel')

                # set values in instance
                for e in org:
                    setattr(djOrg,e,org[e])


                djOrg.clean_Fields()
                validDict = djOrg.validate()
                if validDict:
                    logger.info(f" XX {djOrg} {validDict} ")
                else:
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        logger.info(f" -> {djOrg} ")
                        djOrg.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    else:
                        logger.info(f" >r {djOrg} ")
            else:
                logger.info(f"[Organism] not Found '{org['ORGANISM_NAME']}' ")
                nProc['notFound'] = nProc['notFound'] + 1
        else:
            nProc['notClass'] = nProc['notClass'] + 1
    logger.info(f"[Organism] [{nTotal-(nProc['Saved']+nProc['notClass']+nProc['notFound'])}] {nTotal} {nProc}")

#-----------------------------------------------------------------------------------
def update_Organism_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------
    orgSQL = """
    Select Organism_ID, Organism_Name, Strain_IDS, Strain_Type, Screen_Panel, Strain_Code, Tax_ID, 
        Resistance_Property, Genetic_Property, Growth_Preference, Strain_Origin, Strain_Reference, 
        Strain_Ident, Strain_AddNotes,
        MTA_Document, MTA_Status, Risk_Group, Pathogen, 
        Oxygen_Pref, Lab_Restriction,
        Sequence_Link, Sequence_MLST,
        Source, Source_Code,      
        Biologist
    From Organism
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[Organism] ... ")
    orgLst = OrgDB.get_dict_list(orgSQL)
    nTotal = len(orgLst)
    logger.info(f"[Organism] {nTotal} ")
    OrgDB.close()

    nProc = {}
    nProc['Saved'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0

    # check user
    appuser = ApplicationUser.get(uploaduser)
    # appuser = None
    # if uploaduser:
    #     appuser = get_User_entry(uploaduser)

    for org in orgLst:
        #print(org)
        #print(org['ORGANISM_ID'])
        orgID = reformat_OrganismID(org['ORGANISM_ID'])
        orgClass = orgID.split("_")[0]

        if orgClass in ['GN','GP','FG','MB']:
        # check if instance exists in DB
            djOrg = Organism.get(orgID)
            if djOrg is None:
                djOrg = Organism()
                djOrg.organism_id = orgID

            # check if organism_name exists
            djOrg.organism_name = Taxonomy.get(org['ORGANISM_NAME'],verbose=1)

            if djOrg:
                djOrg.strain_ids = org['STRAIN_IDS']
                djOrg.strain_code= org['STRAIN_CODE']
                djOrg.strain_type = Dictionary.get_DictValues_fromStrList(djOrg.Choice_Dictionary["strain_type"],org['STRAIN_TYPE'],None)
                djOrg.strain_panel = split_StrList(org['SCREEN_PANEL'])
                djOrg.res_property= org['RESISTANCE_PROPERTY']
                djOrg.gen_property= org['GENETIC_PROPERTY']
                djOrg.strain_origin= org['STRAIN_ORIGIN']
                djOrg.reference= org['STRAIN_REFERENCE']
                djOrg.growth_preference= org['GROWTH_PREFERENCE']
                djOrg.strain_notes= org['STRAIN_ADDNOTES']
                djOrg.res_property= org['RESISTANCE_PROPERTY']
                djOrg.source = org['SOURCE']
                djOrg.source_code = org['SOURCE_CODE']
                djOrg.tax_id= org['TAX_ID']
                djOrg.sequence_link= org['SEQUENCE_LINK']
                djOrg.sequence_mlst= org['SEQUENCE_MLST']
                djOrg.strain_identification= org['STRAIN_IDENT']
                djOrg.mta_status = Dictionary.get('License_Status',org['MTA_STATUS'],None)
                djOrg.mta_document= org['MTA_DOCUMENT']
                djOrg.lab_restriction = Dictionary.get(djOrg.Choice_Dictionary["lab_restriction"],org['LAB_RESTRICTION'],None)
                djOrg.risk_group = Dictionary.get(djOrg.Choice_Dictionary["risk_group"],org['RISK_GROUP'],None)
                djOrg.pathogen_group = Dictionary.get(djOrg.Choice_Dictionary["pathogen_group"],None,org['PATHOGEN'])
                djOrg.oxygen_pref = Dictionary.get(djOrg.Choice_Dictionary["oxygen_pref"],org['OXYGEN_PREF'],None)
                djOrg.biologist = ApplicationUser.get(org['BIOLOGIST'])

                djOrg.clean_Fields()
                validDict = djOrg.validate()
                if validDict:
                    logger.info(f" XX {djOrg} {validDict} ")
                else:
                    # --- Upload ---------------------------------------------------------
                    if upload:
                        logger.info(f" -> {djOrg} as {appuser}")
                        djOrg.save(user=appuser)
                        nProc['Saved'] = nProc['Saved'] + 1
                    # else:
                    #     logger.info(f" >r {djOrg} as {appuser}")
            else:
                logger.info(f"[Organism] not Found '{org['ORGANISM_NAME']}' ")
                nProc['notFound'] = nProc['notFound'] + 1
        else:
            nProc['notClass'] = nProc['notClass'] + 1
    logger.info(f"[Organism] [{nTotal-(nProc['Saved']+nProc['notClass']+nProc['notFound'])}] {nTotal} {nProc}")

    


#-----------------------------------------------------------------------------------
def update_Taxonomy_ora(upload=False,uploaduser=None,OutputN=1000):
#-----------------------------------------------------------------------------------
    taxSQL = """
    Select Organism_Name, Organism_Name_Other, Organism_Code, Organism_Class,
        Tax_ID, Parent_Tax_ID, Tax_Rank, Division, Division_Code, Lineage
    From Organism_Tax
    -- Where Organism_Name like 'Klebsiella%'
    """
    OrgDB = oraCastDB.openOrgDB()
    logger.info(f"[Taxonomy] ... ")
    taxLst = OrgDB.get_dict_list(taxSQL)
    nTotal = len(taxLst)
    logger.info(f"[Taxonomy] {nTotal} ")
    OrgDB.close()

    nTime  = zData.Timer(nTotal)
    nProcessed = 0

    checkCharFields = []
    # check user
    appuser = ApplicationUser.get(uploaduser)

    for tax in taxLst:
        
        # set instance
        djTax = Taxonomy.get(tax['ORGANISM_NAME'])
        if djTax is None:
            djTax = Taxonomy()
            djTax.organism_name = tax['ORGANISM_NAME'].strip()
        djTax.other_names = tax['ORGANISM_NAME_OTHER']
        djTax.code = tax['ORGANISM_CODE']
        djTax.org_class = Dictionary.get(djTax.Choice_Dictionary["org_class"],None,tax['ORGANISM_CLASS'],verbose=1)
        djTax.tax_id = tax['TAX_ID']
        djTax.parent_tax_id = tax['PARENT_TAX_ID']
        djTax.tax_rank = tax['TAX_RANK']
        djTax.division = Dictionary.get(djTax.Choice_Dictionary["division"],tax['DIVISION_CODE'],None,verbose=1)            
        if tax['LINEAGE']:
            djTax.lineage = split_StrList(tax['LINEAGE'])
        #djTax.urlname = slugify(tax['ORGANISM_NAME'],lower=False,allow_unicode=False)
        djTax.urlname = slugify(tax['ORGANISM_NAME'],allow_unicode=False)
        
        djTax.clean_Fields()
        validDict = djTax.validate()
        if validDict:
            logger.info(f" XX {djTax} {validDict} ")
        else:
            #--- Upload ---------------------------------------------------------
            if upload:
                #print(f" -> {djTax} as {appuser}")
                djTax.save(user=appuser)
            else:
                logger.info(f" >r {djTax} as {appuser}")

        nProcessed = nProcessed + 1
        if nProcessed%OutputN == 0:
            eTime,sTime = nTime.remains(nProcessed)
            logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} - {djTax} ")

    eTime,sTime = nTime.remains(nProcessed)
    logger.info(f"[Taxonomy] {nTotal} {sTime} - Finished")
