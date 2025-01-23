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
logName = "Update_Readouts"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.StreamHandler()],
    level=logLevel)
logger.info("-------------------------------------------")

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

    
    if prgArgs.plateid :
        tp = TestPlate.get(prgArgs.plateid, WellData=True, verbose=1)
        if tp:
            print(repr(tp))
            _outstr  = f" {tp.plate_quality} Zf: {tp.zfactor:.3f} "
            _outstr += f"[POS: {tp.poscontrol_stats[tp.STATS_MEDIAN]:.3f} {tp.poscontrol_stats[tp.STATS_MAD]:.3f}] "
            _outstr += f"[NEG: {tp.negcontrol_stats[tp.STATS_MEDIAN]:.3f} {tp.negcontrol_stats[tp.STATS_MAD]:.3f}] "
            #_outstr += f"[Edge: {tp.edge_stats[0]:.3f} {tp.edge_stats[1]:.3f}] "
            logger.info(f"[Previous] {tp.plate_id} - {_outstr} ")

            tp.calc_inhibition(verbose=0)

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
