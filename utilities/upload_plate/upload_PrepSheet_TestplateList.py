import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from functools import reduce
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "UploadTestPlateList"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
#    format="[%(name)-20s] %(message)s ",
    format="%(message)s",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def load_PrepSheet_TestPlateList(db,xlPrepSheet,upload=False):
    dfTPl = oraCastDB.read_ExcelSheet(xlPrepSheet,"TestPlateList")
    dfTPl.columns = [c.upper() for c in dfTPl.columns]
    print(dfTPl)

    sCol = {"ASSAY_ID":"ASSAYTYPE_ID",
            "RESULT_TYPE":"RESULT_TYPE",
            "MOTHERPLATE_ID":"MOTHERPLATE_ID",
            "MOTHERPLATE2_ID":"MOTHERPLATE2_ID",
            "PLATING":"PLATING",
            "LABWARE":"LABWARE_ID",
            "TEST_STRAIN":"TEST_STRAIN",
            "TEST_MEDIA":"MEDIA_ID",
            "TEST_DYE":"TEST_DYE",
            "TEST_ADDITIVE":"TEST_ADDITIVE",
            "LAYOUT":"LAYOUT_CONTROL",
            "PROCESSING":"PROCESSING",
            "ISSUES":"ISSUES",
            "SYN_COMPOUNDS_AB":"SYN_COMPOUNDS_A",
            "SYN_COMPOUNDS_POT":"SYN_COMPOUNDS_B"
        }


    if len(dfTPl) > 0:
        for idx, row in dfTPl.iterrows():      
            fUpload = False
            nCnt = oraCastDB.check_PlateID_exists(db,"TestPlate",row['TESTPLATE_ID'])
            if nCnt == 1:
                fUpload = True
            else:
                fUpload = False
                logger.error(f"[CastDB] No TestPlate Found [{row['TESTPLATE_ID']}]")
            if pd.notnull(row['MOTHERPLATE_ID']):
                m1Cnt = oraCastDB.check_PlateID_exists(db,"MasterPlate",row['MOTHERPLATE_ID'])
                if m1Cnt == 1:
                    fUpload = True
                else:
                    fUpload = False
                    logger.error(f"[CastDB] No MotherPlate Found [{row['MOTHERPLATE_ID']}]")
            if pd.notnull(row['MOTHERPLATE2_ID']):
                if row['MOTHERPLATE2_ID']:
                    m2Cnt = oraCastDB.check_PlateID_exists(db,"MasterPlate",row['MOTHERPLATE2_ID'])
                    if m2Cnt == 1:
                        fUpload = True
                    else:
                        fUpload = False
                        logger.error(f"[CastDB] No MotherPlate2 Found [dlTPl[p]['MOTHERPLATE_ID']]")


            if fUpload and upload:
                cTable = "TestPlate"
                sWhere = f" Plate_ID = '{row['TESTPLATE_ID']}' "
                sDict = set_dictFields(row,sCol)
                sSql = db.gen_UpdateSQL(cTable,sDict,sWhere,bindvars=True)

                logger.info(f"[CastDB] Updating {row['TESTPLATE_ID']} PlateData ")
                db.exec(sSql,sDict,commit=True)
            else:
                logger.info(f"[CastDB] NO Updating for {row['TESTPLATE_ID']} ")

    iCol = "TestPlate_ID"
#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from applib.plate.multimode_reader import multimodereader_xls
    from dscreen.models import Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # TestPlate XLSX -------------------------------------------------------------
    if prgArgs.table == 'TestPlateList':



#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [TestPlate]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
    prgParser.add_argument("-e","--excel",default=None,required=True, dest="excelfile", action='store', help="Excel File")
    prgParser.add_argument("--prefix",default=None,required=False, dest="prefix", action='store', help="Prefix to add to PlateID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    from zDjango.djUtils import init_django_dir

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)

#==============================================================================