#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from zUtils import zData

import django
#from djCOADD import djOrgDB
from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

logTime= datetime.datetime.now()
logName = "UploadOrgDB"
logFileName = f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log"

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
    prgParser.add_argument("-f","--vitekfolder",default=None,required=False, dest="vitekfolder", action='store', help="Vitek Folder to parse")
    prgParser.add_argument("-p","--vitekfile",default=None,required=False, dest="vitekfile", action='store', help="Vitek File to parse")
    prgParser.add_argument("--orgbatch",default=None,required=False, dest="orgbatch", action='store', help="OrganismBatch ID")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD"
    uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/Data"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/Data"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/Data"

    xlFiles = {
        'Application': "ApplicationData_v03.xlsx",
        'Drug': "DrugData_v02.xlsx",
    }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    import a_upload_AppUtil as appUtil
    import b_upload_dOrganism as dOrg
    import c_upload_dDrug as dDrug
    import d_upload_Vitek as dVitek

    # Logger ----------------------------------------------------------------
    logger = logging.getLogger(__name__)
    logging.basicConfig(
#    format="%(asctime)s [%(levelname)-8s] [%(name)s] %(message)s ",
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),
              logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logging.INFO)
#    level=logging.DEBUG)

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # Excel  -------------------------------------------------------------
    #uploadFile = os.path.join(uploadDir,"ApplicationData_2022_11_21_JZuegg_v01.xlsx")

    # Table -------------------------------------------------------------
    logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
    logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

    if prgArgs.table == 'User':
        uploadFile = os.path.join(uploadDir,xlFiles['Application'])
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from [User] in {uploadFile}") 
        appUtil.update_AppUser_xls(uploadFile,XlsSheet="User", upload=prgArgs.upload)
    elif prgArgs.table == 'Dictionary':
        uploadFile = os.path.join(uploadDir,xlFiles['Application'])
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Dictionary] in {uploadFile}")
        appUtil.update_Dictionary_xls(uploadFile,XlsSheet="Dictionary", upload=prgArgs.upload,uploaduser=prgArgs.appuser,lower=True) 

    elif prgArgs.table == 'Taxonomy':
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
        dOrg.update_Taxonomy_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    elif prgArgs.table == 'Organism':
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
        dOrg.update_Organism_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    elif prgArgs.table == 'OrgBatch':
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
        dOrg.update_OrgBatch_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    elif prgArgs.table == 'OrgBatchStock':
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
        dOrg.update_OrgBatchStock_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

    elif prgArgs.table == 'Drug':
        uploadFile = os.path.join(uploadDir,xlFiles['Drug'])
        logger.info(f"[Upd_djCOADD] {prgArgs.table} from [Drug] in {uploadFile}") 
        dDrug.update_Drug_xls(uploadFile,XlsSheet="Drug", upload=prgArgs.upload,uploaduser=prgArgs.appuser,lower=True)
    elif prgArgs.table == 'Vitek':
        if prgArgs.vitekfolder:
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from folder {prgArgs.vitekfolder}") 
            dVitek.update_VitekCards(VitekFolder=prgArgs.vitekfolder,upload=prgArgs.upload,uploaduser=prgArgs.appuser)
        if prgArgs.vitekfile:
            logger.info(f"[Upd_djCOADD] {prgArgs.table} from folder {prgArgs.vitekfile}") 
            dVitek.update_VitekCard_single(VitekFile=prgArgs.vitekfile,upload=prgArgs.upload,uploaduser=prgArgs.appuser,OrgBatchID=prgArgs.orgbatch)
    #elif prgArgs.table == 'VitekID':
    #     logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
    #     dDrug.update_VitekID_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    # elif prgArgs.table == 'VitekAST':
    #     logger.info(f"[Upd_djCOADD] {prgArgs.table} from oraCastDB") 
    #     dDrug.update_VitekAST_ora(upload=prgArgs.upload,uploaduser=prgArgs.appuser)

        
    #print(prgArgs)uploadFile,

    #

        
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================