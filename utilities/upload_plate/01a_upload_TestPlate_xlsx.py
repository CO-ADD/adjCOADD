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
logName = "UploadTestPlateXLSX"
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

def main(prgArgs,djDir):

    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from applib.plate.multimode_reader import multimodereader_xls
    from dscreen.models import Screen_Run
    from dsummary.models import Summary_ScreenRun 

    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # TestPlate XLSX -------------------------------------------------------------
    if prgArgs.table == 'TestPlate':
    
        if prgArgs.runid and prgArgs.excelfile:
            new_runid = False
            n_uploads = 0

            djRun = Screen_Run.get(prgArgs.runid)
            if djRun is None:
                djRun = Screen_Run()
                djRun.run_id = prgArgs.runid
                new_runid = True
            
            if new_runid and prgArgs.upload:
                djRun.save()    

            if os.path.isfile(prgArgs.excelfile):
                logger.info(f"[Reading XLSX: {prgArgs.excelfile} ({prgArgs.runid}) ")
                lstTP = multimodereader_xls(prgArgs.excelfile,prgArgs.prefix, verbose=1)

                if prgArgs.upload:
                    _desc = 'TestPlates Saving'
                else:
                    _desc = 'TestPlates Validating'

                for tpDict in tqdm(lstTP, desc=_desc):
                    validStatus = True
                    validDict = {}
                    tpDict['plate'].run_id = djRun
                    
                    tpDict['plate'].set_defaults_model()
                    validDict = tpDict['plate'].validate_model(WellData=False, verbose = 0)
                    if validDict:
                        validStatus = False
                        validDF = pd.DataFrame(validDict)
                        for c in validDF.columns:
                            print(validDF[c].unique())
                        
                    if prgArgs.upload and validStatus:
                        if tpDict['new'] or prgArgs.overwrite:
                            tpDict['plate'].save(verbose=0)
                            n_uploads += 1
            if n_uploads> 0:
                Summary_ScreenRun.update(djRun, to_save=True)

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

    prgParser.add_argument("--django",default='Local',required=True, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
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