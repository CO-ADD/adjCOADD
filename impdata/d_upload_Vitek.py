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

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from apputil.utils import validation_log
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST
from ddrug.utils.import_drug import *

import ddrug.utils.vitek as Vitek

#-----------------------------------------------------------------------------------
def update_VitekPDF(PdfFile=None,VitekFolder=None,ProcessedFolder=None,OrgBatchID=None,
                            upload=False,appuser=None):
#-----------------------------------------------------------------------------------
    lCards,lID,lAST = Vitek.process_VitekPDF(VitekFolder,PdfFile,OrgBatchID=OrgBatchID)

    vLog = validation_log.Validation_Log(PdfFile)

    for c in lCards:
        djCard = imp_VitekCard_fromDict(c,vLog)
        if upload:
            if djCard.VALID_STATUS:
                logger.debug(f" {djCard} {appuser}")
                djCard.save(user=appuser)
            else:
                vLog.add_log('Error','Vitek Card not validated',f"{c['CARD_BARCODE']}",'-')    

    for c in lID:
        djID = imp_VitekID_fromDict(c,vLog)
        if upload:
            if djID.VALID_STATUS:
                logger.debug(f" {djID} {appuser}")
                djID.save(user=appuser)
            else:
                vLog.add_log('Error','Vitek ID not validated',f"{c['CARD_BARCODE']}",'-')    

    for c in lAST:
        djAST = imp_VitekAST_fromDict(c,vLog)
        if upload:
            if djAST.VALID_STATUS:
                logger.debug(f" {djAST} {appuser}")
                djAST.save(user=appuser)
            else:
                vLog.add_log('Error','Vitek AST not validated',f"{c['CARD_BARCODE']} {c['DRUG_NAME']}",'-')    

    if vLog.nLogs['Error'] > 0:
        if upload:
            logger.info(f"[Vitek    ] NOT UPLOADED - {PdfFile}  ")
        vLog.select_unique()
        vLog.show(logTypes=['Error'])
    else:
        if upload:
            logger.info(f"[Vitek    ] -> {ProcessedFolder}  ")
            os.rename(os.path.join(VitekFolder,PdfFile),os.path.join(ProcessedFolder,f"{PdfFile}"))
            #add(cls, LogCode, LogProc,LogType,LogUser,LogObject,LogDesc,LogStatus)
            logDesc = f'{VitekFolder} {PdfFile} [{len(lCards)} {len(lID)} {len(lAST)}]'
            ApplicationLog.add('Import','upload_Vitek','Info',appuser,'Vitek',logDesc,'Completed')

#-----------------------------------------------------------------------------------
def update_VitekCards(VitekFolder=None,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------

    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    # Get PDF Files
    nPdf = 0
    if VitekFolder:
        DirFiles = os.listdir(VitekFolder)
        PdfFiles = [f for f in DirFiles if f.endswith(".pdf")]
        inFiles = [os.path.join(VitekFolder,f) for f in DirFiles if f.endswith(".pdf")]
        nPdf = len(PdfFiles)

    if nPdf>0:
        vCards = []
        vID = []
        vAST = []
        ProcessedFolder = os.path.join(VitekFolder,".Uploaded")
        if not os.path.exists(ProcessedFolder):
            os.makedirs(ProcessedFolder)

        for i in range(nPdf):
            logger.info("-------------------------------------------------------------------------")
            logger.info(f"[Vitek    ] {i+1:3d}/{nPdf:3d} - {PdfFiles[i]}   [{appuser}] ")
            update_VitekPDF(PdfFile=PdfFiles[i],VitekFolder=VitekFolder,ProcessedFolder=ProcessedFolder,OrgBatchID=None,
                            upload=upload,appuser=appuser)
    else:
        logger.info(f"[Vitek    ] NO PDF to process in {VitekFolder}  ")

#-----------------------------------------------------------------------------------
def update_VitekCard_single(VitekFile=None,upload=False,uploaduser=None,OrgBatchID=None):
#-----------------------------------------------------------------------------------

    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    # Get PDF Files
    nPdf = 1
    VitekFolder = os.path.dirname(VitekFile)
    PdfFile = os.path.basename(VitekFile)

    if os.path.exists(VitekFile):
        ProcessedFolder = os.path.join(VitekFolder,".Uploaded")
        if not os.path.exists(ProcessedFolder):
            os.makedirs(ProcessedFolder)

        logger.info("-------------------------------------------------------------------------")
        logger.info(f"[Vitek    ] {PdfFile}   [{appuser}] ")
        update_VitekPDF(PdfFile=PdfFile,VitekFolder=VitekFolder,ProcessedFolder=ProcessedFolder,OrgBatchID=OrgBatchID,
                        upload=upload,appuser=appuser)
    else:
        logger.info(f"[Vitek    ] NO PDF to process in {VitekFolder}  ")

#-----------------------------------------------------------------------------------
def update_xVitekCards(VitekFolder=None,upload=False,uploaduser=None):
#-----------------------------------------------------------------------------------

    appuser = None
    if uploaduser:
        appuser = ApplicationUser.get(uploaduser)

    # Get PDF Files
    nPdf = 0
    if VitekFolder:
        DirFiles = os.listdir(VitekFolder)
        PdfFiles = [f for f in DirFiles if f.endswith(".pdf")]
        inFiles = [os.path.join(VitekFolder,f) for f in DirFiles if f.endswith(".pdf")]
        nPdf = len(PdfFiles)

    if nPdf>0:
        ProcessedFolder = os.path.join(VitekFolder,".Uploaded")
        if not os.path.exists(ProcessedFolder):
            os.makedirs(ProcessedFolder)

        for i in range(nPdf):
            logger.info("-------------------------------------------------------------------------")
            logger.info(f"[Vitek    ] {i+1:3d}/{nPdf:3d} - {PdfFiles[i]}   [{appuser}] ")
            lCards,lID,lAST = Vitek.process_VitekPDF(VitekFolder,PdfFiles[i])

            vLog = Validation_Log(PdfFiles[i])

            for c in lCards:
                djCard = VITEK_Card.check_from_dict(c,vLog)
                if upload:
                    if djCard.VALID_STATUS:
                        logger.debug(f" {djCard} {appuser}")
                        djCard.save(user=appuser)
                    else:
                        vLog.add_log('Error','Vitek Card not validated',f"{c['CARD_BARCODE']}",'-')    

            for c in lID:
                djID = VITEK_ID.check_from_dict(c,vLog)
                if upload:
                    if djID.VALID_STATUS:
                        logger.debug(f" {djID} {appuser}")
                        djID.save(user=appuser)
                    else:
                        vLog.add_log('Error','Vitek ID not validated',f"{c['CARD_BARCODE']}",'-')    

            for c in lAST:
                djAST = VITEK_AST.check_from_dict(c,vLog)
                if upload:
                    if djAST.VALID_STATUS:
                        logger.debug(f" {djAST} {appuser}")
                        djAST.save(user=appuser)
                    else:
                        vLog.add_log('Error','Vitek AST not validated',f"{c['CARD_BARCODE']} {c['DRUG_NAME']}",'-')    

            if vLog.nLogs['Error'] > 0:
                if upload:
                    logger.info(f"[Vitek    ] NOT UPLOADED - {PdfFiles[i]}  ")
                vLog.select_unique()
                vLog.show(logTypes=['Error'])
            else:
                if upload:
                    logger.info(f"[Vitek    ] -> {ProcessedFolder}  ")
                    os.rename(os.path.join(VitekFolder,PdfFiles[i]),os.path.join(ProcessedFolder,f"{PdfFiles[i]}"))
                    #add(cls, LogCode, LogProc,LogType,LogUser,LogObject,LogDesc,LogStatus)
                    logDesc = f'{VitekFolder} {PdfFiles[i]} [{len(lCards)} {len(lID)} {len(lAST)}]'
                    ApplicationLog.add('Import','upload_Vitek','Info',appuser,'Vitek',logDesc,'Completed')

    else:
        logger.info(f"[Vitek    ] NO PDF to process in {VitekFolder}  ")

