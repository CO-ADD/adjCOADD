#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ConvertID"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------

def get_oraProject():
    from oraCastDB.oraCastDB import openCastDB

    fix_StockType = {
        "Master"    : "MasterStock",
        "Master_LN2":"MasterLN2",
        "Stock"     :"Stock",
        "HighUse"   :"HighUse",
    }

    firstStockDate = datetime.date(2010, 1, 1)
    prjSQL = """
     Select PROJECT_ID, PROJECT_NAME, PROJECT_STATUS, PROJECT_TYPE, LIBRARY_NAME,


        COADD_ID, CPOZ_ID,
        RECEIVED, COMPLETED,
        COMPOUND_COMMENT, COMPOUND_STATUS,
        GROUP_ID, CONTACT_A_ID, CONTACT_B_ID,
        ORGANISATION, COUNTRY,
        REPORT_COMMENT, REPORT_HC_DATE, REPORT_HV_DATE, REPORT_PS_DATE, REPORT_STATUS,
        SCREEN_COMMENT,
        SCREEN_CONC, SCREEN_CONC_UNIT, SCREEN_STATUS,
        STOCK_COMMENT, STOCK_CONC, STOCK_CONC_UNIT, STOCK_CONTAINER, STOCK_STATUS,
        DATA_COMMENT, DATA_STATUS,
        ANTIMICRO_STATUS,
        PROJECT_ACTION, PROJECT_COMMENT,
        PROVIDED_COMMENT, PROVIDED_CONTAINER,
        PUB_DATE, PUB_STATUS,
    From Project
    -- Where Organism_Name like 'Klebsiella%'
    """
    CastDB = openCastDB()
    logger.info(f"[Projects] ... ")
    prjList = CastDB.get_dict_list(prjSQL)
    nTotal = len(prjList)
    logger.info(f"[Projects] {nTotal} ")
    CastDB.close(prjList)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dsample.models import Project
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------
    runTables = ["CompoundID","ProjectID"]

    if prgArgs.table == "ProjectID" :

        prjList = get_oraProject()




        # ExcelFile = prgArgs.file
        # SheetName = "ProjectID"
        # print(f"[Reading Excel] {ExcelFile} {SheetName} ")
        # prjDF = pd.read_excel(ExcelFile,sheet_name=SheetName)
        # prjDF.project_name = prjDF.project_name.fillna('')
        
        # for idx,row in tqdm(prjDF.iterrows(), total=prjDF.shape[0]):
        #     djE = Convert_ProjectID.get(row['ora_project_id'])
        #     if djE is None:
        #         djE = Convert_ProjectID()
        #         djE.ora_project_id = row['ora_project_id']
        #         djE.project_id = row['project_id']
        #         djE.project_name = row['project_name']
        #         if prgArgs.upload:
        #             djE.clean_Fields()
        #             djE.save()
            # else:
            #     print(f"[Exists already] {row['ora_project_id']} {row['project_id']} ")

    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.config == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.config == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.config == 'Laptop':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
