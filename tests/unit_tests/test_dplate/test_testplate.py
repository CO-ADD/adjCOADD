import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from functools import reduce
from pathlib import Path

from tqdm import tqdm

from zSql import zSqlConnector
from zDjango.djUtils import init_django_dir
import django

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ActStrDR"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.StreamHandler()],
    level=logLevel)
logger.info("-------------------------------------------")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):
    print(djDir['Project'],djDir['djPrj'])
    sys.path.append(djDir['djPrj'])
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", f"{djDir['Project']}.settings")
    django.setup()

    logging.getLogger().addHandler(logging.FileHandler(logFileName,mode='w'))
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Data    : {djDir['dataDir']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    from dplate.models import TestPlate

    existingPlateID = 'TP00360-04C'
    print("")
    print("------------------------------------------------")
    print(f"-- Test Existing TestPlate {existingPlateID}")
    tp = TestPlate.get(existingPlateID, WellData=True, verbose=1)
    print(repr(tp))
    # tp.calc_inhibition(verbose=1)
    # print(f"[PosCtrl] {tp.poscontrol_stats}")
    # print(f"[NegCtrl] {tp.negcontrol_stats}"),
    # print(f"[Sample ] {tp.sample_stats}")
    # print(f"[Edge   ] {tp.edge_stats}")
    # print(f"[ZFactor] {tp.zfactor}")
    
    # for w in tp.wells:
    #     print(f"{tp.wells[w]} {tp.wells[w].inhibition}")

    #--------------------------------------------------------
    newPlateID = 'xx_test'
    print("")
    print("------------------------------------------------")
    print(f"-- Test TestPlate.new {newPlateID}")
    np = TestPlate.new('xx_test',96)
    if np is None:
        np = TestPlate.get(newPlateID.upper(), WellData=True, verbose=1)
    
    if np:
        print(repr(np))
        # for w in np.wells:
        #     print(f"{np.wells[w]} {np.wells[w].inhibition}")

        #print(f" (get_welldata) {np.get_welldata()}")
        np.set_defaults_model()
        np.calc_inhibition(verbose=1)
        if prgArgs.upload:
            np.save()

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")

    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=False, dest="table", action='store', help="Table to upload [CompoundID]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=True, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

