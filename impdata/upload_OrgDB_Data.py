#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

# from zUtils import zData

import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

import logging
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_OrgDB_Data', 
                                description="Uploading data to adjCOADD from Oracle or Excel")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--orgbatch",default=None,required=False, dest="orgbatch", action='store', help="OrganismBatch ID")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD"
    uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/impdata/Data"
    orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/impdata/Data"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/impdata/Data"

    xlFiles = {
        'Application': "ApplicationData_v05.xlsx",
        'Drug': "DrugData_v04.xlsx",
        'MIC': "LMIC_Data_v06.xlsx",
        'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    import a_upload_AppUtil as appUtil
    import b_upload_dOrganism as dOrg
    import c_upload_dDrug as dDrug
    import d_upload_Vitek as dVitek
    import e_upload_MIC as dMIC
    import f_upload_gene as dGene

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "UploadOrgDB"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
        level=logLevel)
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # Excel  -------------------------------------------------------------
    #uploadFile = os.path.join(uploadDir,"ApplicationData_2022_11_21_JZuegg_v01.xlsx")

    # Table -------------------------------------------------------------

    choiceTables = ['User','Dictionary',
                    'Taxonomy','Organism','OrgBatch','OrgBatchStock','OrgBatchImages','OrgCulture',
                    'Drug','MICPub','MICCollab','MICCOADD','BP',
                    'Vitek',
                    'wgsAssembly','wgsFastA','wgsAMR',
                    ]
    if prgArgs.table in choiceTables:

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

        if prgArgs.table == 'User':
            prgArgs.upload = False
            uploadFile = os.path.join(uploadDir,xlFiles['Application'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [User] in {uploadFile} [Upload: {prgArgs.upload}]") 
            appUtil.update_AppUser_xls(uploadFile,XlsSheet="User", upload=prgArgs.upload)

        elif prgArgs.table == 'Dictionary':
            prgArgs.upload = False
            uploadFile = os.path.join(uploadDir,xlFiles['Application'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Dictionary] in {uploadFile} [Upload: {prgArgs.upload}]")
            appUtil.update_Dictionary_xls(uploadFile,XlsSheet="Dictionary", upload=prgArgs.upload,uploaduser=prgArgs.appuser,lower=True) 

        elif prgArgs.table == 'Taxonomy':
            prgArgs.upload = False
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraOrgDB [Upload: {prgArgs.upload}]") 
            dOrg.update_Taxonomy_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'Organism':
            prgArgs.upload = False
            uploadFile = os.path.join(orgdbDir,xlFiles['OrgDB'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Organism] in {uploadFile} [Upload: {prgArgs.upload}]") 
            dOrg.update_Organism_xls(uploadFile,XlsSheet="Organism",upload=prgArgs.upload,uploaduser=prgArgs.appuser)

            #logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraOrgDB") 
            #dOrg.update_Organism_ora(uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'OrgBatch':
            prgArgs.upload = False
            uploadFile = os.path.join(orgdbDir,xlFiles['OrgDB'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [OrgBatch] in {uploadFile} [Upload: {prgArgs.upload}]") 
            dOrg.update_OrgBatch_xls(uploadFile,XlsSheet="OrgBatch",upload=prgArgs.upload,uploaduser=prgArgs.appuser)

            #logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraOrgDB") 
            #dOrg.update_OrgBatch_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'OrgBatchStock':
            prgArgs.upload = False
            uploadFile = os.path.join(orgdbDir,xlFiles['OrgDB'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [OrgBatch Stock] in {uploadFile} - Upload: {prgArgs.upload}")
            dOrg.delete_OrgBatchStock(upload=prgArgs.upload)
            dOrg.update_OrgBatchStock_xls(uploadFile,XlsSheet="Stock",upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'OrgBatchImages':
            #prgArgs.upload = False
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from folder {prgArgs.directory} [Upload: {prgArgs.upload}]") 
            dOrg.update_OrgBatchImg(prgArgs.directory,upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'OrgCulture':
            prgArgs.upload = False
            uploadFile = os.path.join(orgdbDir,xlFiles['OrgDB'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Organims Culture] in {uploadFile} [Upload: {prgArgs.upload}]")
            dOrg.delete_OrgCulture(upload=prgArgs.upload)
            dOrg.update_OrgCulture_xls(uploadFile,XlsSheet="Culture",upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'Drug':
            prgArgs.upload = False
            uploadFile = os.path.join(uploadDir,xlFiles['Drug'])
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Drug] in {uploadFile} [Upload: {prgArgs.upload}]") 
            dDrug.update_Drug_xls(uploadFile,XlsSheet="Drug", upload=prgArgs.upload,uploaduser=prgArgs.appuser,lower=True)

        elif prgArgs.table == 'Vitek':
            if prgArgs.directory:
                logger.info(f"[Upd_djCOADD] {prgArgs.table} from folder {prgArgs.directory} [Upload: {prgArgs.upload}]") 
                dVitek.update_VitekCards(VitekFolder=prgArgs.directory,upload=prgArgs.upload,uploaduser=prgArgs.appuser)
            if prgArgs.file:
                logger.info(f"[Upd_djCOADD] {prgArgs.table} from file {prgArgs.file} [Upload: {prgArgs.upload}]") 
                dVitek.update_VitekCard_single(VitekFile=prgArgs.file,upload=prgArgs.upload,uploaduser=prgArgs.appuser,OrgBatchID=prgArgs.orgbatch)

        elif prgArgs.table == 'MICPub':
            prgArgs.upload = False
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraOrgDB [Upload: {prgArgs.upload}]") 
            dMIC.update_MICPub_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'MICCollab':
            prgArgs.upload = False
            uploadFile = os.path.join(uploadDir,xlFiles['MIC']) 
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from [MIC] in {uploadFile} [Upload: {prgArgs.upload}]")
            dMIC.update_MICPub_xls(uploadFile,XlsSheet="MIC",upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'MICCOADD':
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB {prgArgs.runid} [Upload: {prgArgs.upload}]")
            if  prgArgs.runid:
                dMIC.update_MICCOADD_ora(prgArgs.runid,upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'BP':
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraOrgDB [Upload: {prgArgs.upload}]")
            dMIC.update_Breakpoints_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)


        elif prgArgs.table == 'wgsAssembly':
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from 02_Assembly {prgArgs.runid} [Upload: {prgArgs.upload}]")
            #dGene.update_WGSCOADD_Trim(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
            dGene.update_WGSCOADD_Assembly(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
            
        elif prgArgs.table == 'wgsFastA':
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from 03_FastA {prgArgs.runid} [Upload: {prgArgs.upload}]")
            dGene.update_WGSCOADD_FastA(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        elif prgArgs.table == 'wgsAMR':
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from 03_FastA {prgArgs.runid} [Upload: {prgArgs.upload}]")
            #dGene.update_WGSCOADD_AMR(Methods= ['AMR Finder'],upload=prgArgs.upload,uploaduser=prgArgs.appuser)
            dGene.update_WGSCOADD_AMR(Methods= ['Abricate card'],upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    else:
        logger.error(f"[Upd_djCOADD] {prgArgs.table} not in {choiceTables}")
                     
    logger.info(f"-------------------------------------------------------------------")
    logger.info(f"LogFile        : {logFileName}")

    #

        
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
