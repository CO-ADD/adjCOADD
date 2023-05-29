import os,sys
import shutil
import datetime
import csv
import datetime
from dateutil import parser
import pandas as pd

import pdfplumber
import re
from django.core.cache import cache

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from apputil.utils.validation_log import Validation_Log
from ddrug.utils.import_drug import imp_VitekCard_fromDict,imp_VitekAST_fromDict,imp_VitekID_fromDict

import logging
logger = logging.getLogger(__name__)
__version__ = "1.1"

#-----------------------------------------------------------------------------------
def upload_VitekPDF_Process(Request, DirName, FileList, OrgBatchID=None,upload=False,appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single Vitek PDF, given by:
        Request : Objects to pass state through the system, including user model instance: e.g., request.user
        DirName : FolderName
        PdfName : PdfName without FolderName
        OrgBatchID: Optional to overwrite OrgBatchID from PDF
        upload : Validation only (False) or Validation and Upload (True)
        appuser : User Instance of user uploading

    """

    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0

    # Get PDF Files in case none given
    if nFiles == 0:
        if DirName:
            DirFiles = os.listdir(DirName)
            FileList = [f for f in DirFiles if f.endswith(".pdf")]
            nFiles = len(FileList)

    valLog = Validation_Log("upload_VitekPDF_List")

    if nFiles > 0:
        for i in range(nFiles):
            logger.info(f"[upload_VitekPDF_List] {i+1:3d}/{nFiles:3d} - {FileList[i]}   [{appuser}] ")
            upload_VitekPDF(DirName,FileList[i],OrgBatchID=OrgBatchID,upload=upload,appuser=appuser,valLog=valLog)

    else:
        logger.info(f"[upload_VitekPDF_List] NO PDF to process in {DirName}  ")

    valLog.select_unique()
    
    return(valLog)

#-----------------------------------------------------------------------------------
def upload_VitekPDF_List(DirName, FileList, OrgBatchID=None,upload=False,appuser=None):
#-----------------------------------------------------------------------------------
    """
    Uploads (upload=True) the data from a single Vitek PDF, given by:
        DirName : FolderName
        PdfName : PdfName without FolderName
        OrgBatchID: Optional to overwrite OrgBatchID from PDF
        upload : Validation only (False) or Validation and Upload (True)
        appuser : User Instance of user uploading

    """
    if FileList:
        nFiles = len(FileList)
    else:
        nFiles = 0

    # Get PDF Files in case none given
    if nFiles == 0:
        if DirName:
            DirFiles = os.listdir(DirName)
            FileList = [f for f in DirFiles if f.endswith(".pdf")]
            nFiles = len(FileList)

    valLog = Validation_Log("upload_VitekPDF_List")
    if nFiles > 0:
        for i in range(nFiles):
            # Cancel Process, If the cancel flag is set, break the loop, 
            # For break a running upload
            logger.info(f"[upload_VitekPDF_List] {i+1:3d}/{nFiles:3d} - {FileList[i]}   [{appuser}] ")
            upload_VitekPDF(DirName,FileList[i],OrgBatchID=OrgBatchID,upload=upload,appuser=appuser,valLog=valLog)
    else:
        logger.info(f"[upload_VitekPDF_List] NO PDF to process in {DirName}  ")

    valLog.select_unique()    
    return(valLog)

#-----------------------------------------------------------------------------
def upload_VitekPDF(DirName,FileName,OrgBatchID=None,upload=False,appuser=None,valLog=None):
    """
    Uploads (upload=True) the data from a single Vitek PDF, given by:
        DirName : FolderName
        PdfName : PdfName without FolderName
        OrgBatchID: Optional to overwrite OrgBatchID from PDF
        upload : Validation only (False) or Validation and Upload (True)
        appuser : User Instance of user uploading
        volLog: Validation_Log to collect ['Info','Warning','Error']

    """
#-----------------------------------------------------------------------------
    nProc = {}
    nProc['Saved'] = 0
    nProc['notValid'] = 0

    if not valLog:
        valLog = Validation_Log(FileName)

    lCards,lID,lAST = process_VitekPDF(DirName,FileName,OrgBatchID=OrgBatchID)

    nProcess = {'Uploaded':0,'Invalid':0}
    for c in lCards:
        djCard = imp_VitekCard_fromDict(c,valLog,upload=upload)
        if upload:
            if djCard.VALID_STATUS:
                logger.debug(f" {djCard} <- {FileName}")
                djCard.save(user=appuser)
                nProc['Saved'] = nProc['Saved'] + 1
            else:
                # nFailed += 1
                valLog.add_log('Error','Vitek Card not validated',f"{c['CARD_BARCODE']}",'-')    
                nProc['notValid'] = nProc['notValid'] + 1

    for c in lID:
        djID = imp_VitekID_fromDict(c,valLog,upload=upload)
        if upload:
            if djID.VALID_STATUS:
                logger.debug(f" {djID} <- {FileName}")
                djID.save(user=appuser)
                nProc['Saved'] = nProc['Saved'] + 1
            else:
                valLog.add_log('Error','Vitek ID not validated',f"{c['CARD_BARCODE']}",'-')    
                nProc['notValid'] = nProc['notValid'] + 1

    for c in lAST:
        djAST = imp_VitekAST_fromDict(c,valLog,upload=upload)
        if upload:
            if djAST.VALID_STATUS:
                logger.debug(f" {djAST} <- {FileName}")
                djAST.save(user=appuser)
                nProc['Saved'] = nProc['Saved'] + 1
            else:
                valLog.add_log('Error','Vitek AST not validated',f"{c['CARD_BARCODE']} {c['DRUG_NAME']}",'-')    
                nProc['notValid'] = nProc['notValid'] + 1

    # if valLog.nLogs['Error'] > 0:
    #     if upload:
    #         logger.info(f"[upload_VitekPDF] NOT UPLOADED - {FileName}  ")
    #     valLog.select_unique()
    #     valLog.show(logTypes=['Error'])
    else:
        if upload:
            logDesc = f'{FileName} [Cards:{len(lCards)} ID:{len(lID)} AST:{len(lAST)}] : {nProc}'
            ApplicationLog.add('Import','upload_VitekPDF','Info',appuser,'Vitek',logDesc,'Completed')
            # if nProc['notValid'] == 0:
            #   removeFile(DirName,FileName)

    valLog.select_unique()    
    return(valLog)

#-----------------------------------------------------------------------------
def process_VitekPDF(DirName,PdfName,OrgBatchID=None):
    """
    Reads a single Vitek PDF extraxts the information into
    Lists of Cards, AST and ID dictionaries

    """
#-----------------------------------------------------------------------------
    pVitek = parse_VitekPDF(DirName,PdfName,OrgBatchID=None)

    lstCards = []
    lstID = []
    lstAST = []
    for pv in pVitek:
        k = pv.keys()
        #print(f"\n**\n {pv}\n**\n")

        # if any('Barcode' in x for x in k):
        # #if 'ID_Card_Barcode' in pv or 'AST_Card_Barcode' in pv:
        #     lstCards.append(dict_Vitek_Card(pv))
        #     print("Card")

        if 'ID' in k :
            if pv['ID']:
                logger.info(f"[Vitek-ID ]  {pv['OrgBatchID']} - {pv['ID_Card']:10s} ({pv['ID_Card_Barcode']}) - {pv['Organism']} ")
                xCard = dict_Vitek_Card(pv,'ID')
                for x in xCard:
                    lstCards.append(x)
                xID = dict_Vitek_ID(pv)
                for x in xID:
                    lstID.append(x)

        if 'AST' in k :
            if pv['AST']:
                logger.info(f"[Vitek-AST]  {pv['OrgBatchID']} - {pv['AST_Card']:10s} ({pv['AST_Card_Barcode']}) - {pv['Organism']} ")
                xCard = dict_Vitek_Card(pv,'AST')
                for x in xCard:
                    lstCards.append(x)
                xAST = dict_Vitek_AST(pv)
                for x in xAST:
                    lstAST.append(x)
    logger.info(f"[Vitek    ] Cards: {len(lstCards):4d} - ID: {len(lstID):4d} - AST: {len(lstAST):4d} ")
    return(lstCards,lstID,lstAST)

#-----------------------------------------------------------------------------
def parse_VitekPDF(DirName,PdfName,OrgBatchID=None):
#-----------------------------------------------------------------------------
    lstVitek = []
    logging.getLogger().setLevel(logging.WARNING)
    with pdfplumber.open(os.path.join(DirName,PdfName)) as pdf:
        nPage = 0
        df = {}
        for page in pdf.pages:
            #print(df)
            df['DirName'] = DirName
            df['FileName'] = PdfName
            if_IsolateData = False
            nPage += 1

            # - Header --------------------------------------------------------------------------------
            txt_lst = page.extract_text().splitlines()
            for l in txt_lst:
                mX = re.search('Isolate:(.+?)\((.+?)\)', l)
                if mX:
                    if_IsolateData = True
                    xName = mX.group(1)
                    df['VitekProcess'] = mX.group(2)
                    mOrg = re.search('([a-zA-Z]+)(\d+)B(\d+)(.+)',xName)
                    if mOrg:
                        df['OrganismID'] = mOrg.group(1)+f"_{int(mOrg.group(2)):04d}"
                        df['BatchID'] = int(mOrg.group(3))
                        df['OrgBatchID'] = df['OrganismID']+f"_{int(df['BatchID']):02d}"
                    else:
                        xLst = xName.strip().split('-')
                        df['OrganismID'] = xLst[0][0:2]+f"_{int(xLst[0][2:]):04d}"
                        df['BatchID'] = int(xLst[1])
                        df['OrgBatchID'] = df['OrganismID']+f"_{int(df['BatchID']):02d}"
                    if OrgBatchID is not None:
                        # Overwrite OrgBatchID
                        xLst = OrgBatchID.strip().split('_')
                        df['OrganismID'] = "_".join(xLst[0:2])
                        df['BatchID'] = int(xLst[2])
                        df['OrgBatchID'] = df['OrganismID']+f"_{int(df['BatchID']):02d}"
                        
            # If Header contains Isolate/Card data
            if if_IsolateData:
                to_be_Saved = False
                if_VitekID = False
                prev_col1 = ''

                for l in txt_lst:
                    if 'Card Type:' in l:
                        xl = l.split(' ')
                        if 'AST' in xl[2]:
                            df['AST_Card'] = xl[2]
                            df['AST_Card_Barcode'] = xl[5]
                            df['AST_Instrument'] = " ".join(xl[8:])
                        else:    
                            df['ID_Card'] = xl[2]
                            df['ID_Card_Barcode'] = xl[5]
                            df['ID_Instrument'] = " ".join(xl[8:])

                # Iterate through all the tables on that page
                for table in page.extract_tables():
                    if_MIC_Section = False
                    if_ID_Section = False
                    if_AST_Section = False
                    prevRow1 = ""
                    prevRow4 = ""
                    for row in table:
                        #print(row)
                        if row[0]:
                            col1 = row[0].replace("\n", "")
                        else:
                            col1 = prev_col1
                            prev_col1 = ''

                        # - Comment Table -----------------------------
                        #if 'Comments:' == col1:
                        #    to_be_Saved = True

                        # - Indentification Section (2 rows) -----------------------------
                        if if_ID_Section:
                            df['ID_Analysis'] = row[2].replace('Analysis Time: ','')
                            if_ID_Section = False
                        if 'Identification Information' == col1:
                            if row[1] != '':
                                if_ID_Section = True

                                #df['ID_Card'] = row[1].split(' ')[1]
                                df['ID_Card_LotN'] = row[2].split(' ')[1]
                                dt = parser.parse(row[3].replace('Expires: ',''),tzinfos={'AEST':10 * 3600})
                                df['ID_Expiry'] = dt.date().strftime('%Y-%m-%d')
                                prev_col1 = 'ID'

                        # - Organism Origin Section (1 row)-----------------------------
                        if 'Organism Origin' == col1:
                            df['Organism_Origin'] = row[1]
                            if 'VITEK 2' == row[1]:
                                if_VitekID = True

                        # - Selected Organism Section (2 rows) -----------------------------
                        if 'Selected Organism' == col1:
                            if if_VitekID:
                                sOrg = {}
                                xcell = [x.split(' ') for x in row[1].splitlines()]
                                #print(xcell)
                                if len(xcell[0])>2:
                                    if 'Probability' in xcell[0][1]:
                                        sOrg['Organism'] = " ".join(xcell[0][2:])
                                        sOrg['Probability'] = xcell[0][0]
                                    else:
                                        sOrg['Organism'] = " ".join(xcell[0])
                                        sOrg['Probability'] = 0

                                else:
                                    sOrg['Organism'] = xcell[0][0] + " " + xcell[0][1]
                                    sOrg['Probability'] = 0
                                if len(xcell[1])>2:
                                    if 'Confidence' in xcell[1][2]:
                                        sOrg['Confidence'] = " ".join(xcell[1][3:])
                                    else:
                                        sOrg['Confidence'] = "-"

                                #print(sOrg['Organism'])
                                df['ID'] = sOrg
                                df['Organism'] = sOrg['Organism']
                            else:
                                df['ID'] = None
                                df['Organism'] = row[1]

                       # - Analysis Organisms and Tests to Separate Section (n rows) -----------------------------
                        #print(f"{row} -> {col1}")
                        if 'Contraindicating Typical Biopattern' in col1:
#                        if 'Analysis Organisms and Tests to Separate' in col1:
                            lowLst =[]
                            for r in row[0].split('\n'):
                                lowLst.append(" ".join(r.split(',')[0].split(' ')[:2]))
                            if len(lowLst)>1:
                                df['ID']['Organism'] = ", ".join(lowLst[1:])

                        # - Susceptibility Information Section (2 rows) -----------------------------
                        if if_AST_Section:
                            df['AST_Analysis'] = row[2].replace('Analysis Time: ','')
                            if_AST_Section = False     
                        if 'Susceptibility Information' == col1:
                            #to_be_Saved = True
                            if_AST_Section = True 
                            #df['AST_Card'] = row[1].split(' ')[1]
                            df['AST_Card_LotN'] = row[2].split(' ')[1]
                            df['AST'] = {}
                            dt = parser.parse(row[3].replace('Expires: ',''),tzinfos={'AEST':10 * 3600})
                            df['AST_Expiry'] = dt.date().strftime('%Y-%m-%d')
                            prev_col1 = 'AST'

                        if 'AST' == col1:
                            dt = parser.parse(row[3].replace('Completed: ',''),tzinfos={'AEST':10 * 3600})
                            df['AST_Date'] = dt.date().strftime('%Y-%m-%d')
                        if 'ID' == col1:
                            dt = parser.parse(row[3].replace('Completed: ',''),tzinfos={'AEST':10 * 3600})
                            df['ID_Date'] = dt.date().strftime('%Y-%m-%d')

                        # - AES Findings Section --------------------------------------------    
                        if 'AES Findings:' == col1:
                            if 'EUCAST' in row[1]:
                                df['AST_Interpretation'] = 'EUCAST'
                            elif 'CLSI' in row[1]:
                                df['AST_Interpretation'] = 'CLSI'
                            else:
                                df['AST_Interpretation'] = 'Other'       
                            if_MIC_Section = False
                            to_be_Saved = True

                        # - MIC Section -----------------------------------------------------   
                        if if_MIC_Section:                            
                            if row[1] != '':
                                if row[0].replace("\n", "") in ['Urine','Other']:
                                    # fix multiple entries for one drug [Urine Others] taking the first instance
                                    if prevRow1 != "":
                                        df['AST'][prevRow1] = {'MIC': row[1], 'Interpretation': row[2]}
                                else:
                                    df['AST'][row[0].replace("\n", "")] = {'MIC': row[1], 'Interpretation': row[2]}
                                prevRow1 = ""
                            else:
                                prevRow1 = row[0].replace("\n", "")

                            if row[4] != '':
                                if row[3].replace("\n", "") in ['Urine','Other']:
                                    # fix multiple entries for one drug [Urine Others] taking the first instance
                                    if prevRow4 != "":
                                        df['AST'][prevRow4] = {'MIC': row[4], 'Interpretation': row[5]}
                                else:
                                    df['AST'][row[3].replace("\n", "")] = {'MIC': row[4], 'Interpretation': row[5]}
                                prevRow4 = ""
                            else:
                                prevRow4 = row[0].replace("\n", "")

                        if 'Antimicrobial' == col1:
                            if_MIC_Section = True

                #print(f"{to_be_Saved} - {df}")
                if to_be_Saved:
                    df['PageNo'] = nPage
                    lstVitek.append(df)
                    df = {}
            else:
                # Last Isolate page and still Data
                if df:
                    df['PageNo'] = nPage - 1
                    lstVitek.append(df)
                    df = {}

        # Last Page and still Isolate Data
        if if_IsolateData:
            if df:
                df['PageNo'] = nPage - 1
                lstVitek.append(df)
                df = {}

    logging.getLogger().setLevel(logging.INFO)
    #print(lstVitek) 
    return(lstVitek)

# --------------------------------------------------------------------
def dict_Vitek_Card(pCard,card_type):
# --------------------------------------------------------------------
    if card_type == 'AST':
        if 'AST_Card_Barcode' in pCard:
            dfCard = {}
            dfCard['CARD_TYPE'] = 'AST'
            dfCard['CARD_CODE'] = pCard['AST_Card']
            dfCard['CARD_BARCODE'] = pCard['AST_Card_Barcode']
            dfCard['INSTRUMENT'] = pCard['AST_Instrument']
            dfCard['FILENAME'] = pCard['FileName']
            dfCard['EXPIRY_DATE'] = pCard['AST_Expiry']
            dfCard['PROCESSING_DATE'] = pCard['AST_Date']
            dfCard['ANALYSIS_TIME'] = pCard['AST_Analysis']
            dfCard['ORGANISM_ID'] = pCard['OrganismID']
            dfCard['BATCH_ID'] = pCard['BatchID']
            dfCard['ORGBATCH_ID'] = pCard['OrgBatchID']
            return([dfCard])         
        else:
            return([])
    if card_type == 'ID':
        if 'ID_Card_Barcode' in pCard:
            dfCard = {}
            dfCard['ORGANISM_ID'] = pCard['OrganismID']
            dfCard['BATCH_ID'] = pCard['BatchID']
            dfCard['ORGBATCH_ID'] = pCard['OrgBatchID']         
            dfCard['CARD_TYPE'] = 'ID'
            dfCard['CARD_CODE'] = pCard['ID_Card']
            dfCard['CARD_BARCODE'] = pCard['ID_Card_Barcode']
            dfCard['INSTRUMENT'] = pCard['ID_Instrument']
            dfCard['FILENAME'] = pCard['FileName']
            if 'ID_Expiry' in pCard:
                dfCard['EXPIRY_DATE'] = pCard['ID_Expiry']
                dfCard['PROCESSING_DATE'] = pCard['ID_Date']
                dfCard['ANALYSIS_TIME'] = pCard['ID_Analysis']
            return([dfCard])         
        else:
            return([])
    return([])

# --------------------------------------------------------------------
def dict_Vitek_ID(pCard):
# --------------------------------------------------------------------
    dfID = {}
    dfID['CARD_CODE'] = pCard['ID_Card']
    dfID['CARD_BARCODE'] = pCard['ID_Card_Barcode']
    dfID['VITEK_PROCESS'] = pCard['VitekProcess']
    dfID['FILENAME'] = pCard['FileName']
    dfID['PAGENO'] = pCard['PageNo']
    dfID['ORGANISM_ID'] = pCard['OrganismID']
    dfID['BATCH_ID'] = pCard['BatchID']
    dfID['ORGBATCH_ID'] = pCard['OrgBatchID']         
    xID = pCard['ID']
    dfID['ID_ORGANISM'] = xID['Organism']
    dfID['ID_PROBABILITY'] = xID['Probability']
    dfID['ID_CONFIDENCE'] = xID['Confidence'].replace(' identification','')
    return([dfID])


# --------------------------------------------------------------------
def dict_Vitek_AST(pCard):
# --------------------------------------------------------------------
    xAST = pCard['AST']
    lAST = []
    for pAST in xAST:
        dfAST = {}
        dfAST['CARD_CODE'] = pCard['AST_Card']
        dfAST['CARD_BARCODE'] = pCard['AST_Card_Barcode']
        dfAST['VITEK_PROCESS'] = pCard['VitekProcess']
        dfAST['FILENAME'] = pCard['FileName']
        dfAST['PAGENO'] = pCard['PageNo']
        dfAST['ORGANISM_ID'] = pCard['OrganismID']
        dfAST['SELECTED_ORGANISM'] = pCard['Organism']
        dfAST['ORGANISM_ORIGIN'] = pCard['Organism_Origin']
        dfAST['BATCH_ID'] = pCard['BatchID']
        dfAST['ORGBATCH_ID'] = pCard['OrgBatchID']         

        xComment = []
        #dfAST['BP_COMMENT'] = ''
        dfAST['BP_PROFILE'] = xAST[pAST]['Interpretation']
        dfAST['BP_SOURCE'] = pCard['AST_Interpretation']

        nDrug = pAST.replace('/','|').replace(' Acid',' acid')
        if '#' in nDrug:
            nDrug = nDrug.replace('#','')
            xComment.append('Disabled bioART Limitation Rule')
            #dfAST['BP_COMMENT'] = 'Disabled bioART Limitation Rule'
        dfAST['DRUG_NAME'] = nDrug

        dfAST['MIC'] = xAST[pAST]['MIC']
        if dfAST['MIC'] == 'TRM':
            xComment.append('Insufficient incubation time for analysis')
            #dfAST['BP_COMMENT'] = 'Insufficient incubation time for analysis'
            dfAST['MIC'] = ''
        elif dfAST['MIC'] == '(-)':
            xComment.append('Susceptibility testing not recommended')
            #dfAST['BP_COMMENT'] = 'Susceptibility testing not recommended'
            dfAST['MIC'] = ''
            dfAST['BP_PROFILE'] = ''
        elif '**' in dfAST['MIC']:
            xComment.append('User modified')
            dfAST['BP_COMMENT'] = 'User modified'
            #dfAST['MIC'] = dfAST['MIC'].replace('**','')
        elif '*' in dfAST['MIC']:
            xComment.append('AES modified')
            #dfAST['BP_COMMENT'] = 'AES modified'
            dfAST['MIC'] = dfAST['MIC'].replace('*','')

        if dfAST['DRUG_NAME'] == 'ESBL':
            xComment.append('(FEP 1, CTX 0.5, CAZ 0.5, FEP/CA 1/10, CTX/CA 0.5/4, CAZ/CA 0.5/4)')
            #dfAST['BP_COMMENT'] = dfAST['BP_COMMENT'] + " (FEP 1, CTX 0.5, CAZ 0.5, FEP/CA 1/10, CTX/CA 0.5/4, CAZ/CA 0.5/4)"
        elif dfAST['DRUG_NAME'] == 'Inducible Clindamycin Resistance':
            xComment.append('(CM 0.5, CM/E 0.25/0.5)')
            #dfAST['BP_COMMENT'] = dfAST['BP_COMMENT']+ " (CM 0.5, CM/E 0.25/0.5)"

        if len(xComment)>0:
            dfAST['BP_COMMENT'] = "; ".join(xComment)
        else:
            dfAST['BP_COMMENT'] = ''
        #print(dfAST)
        lAST.append(dfAST)
    return(lAST)



