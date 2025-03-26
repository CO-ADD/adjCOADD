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
logName = "Update_BPProfile"
logDir = "log"
logFileName = os.path.join(logDir,f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

if not os.path.isdir(logDir):
    os.mkdir(logDir)

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
def main(prgArgs,djDir):
#-----------------------------------------------------------------------------------

    django.setup()

    # from apputil.models import Dictionary
    # from dsample.models import COADD_Compound, Compound_Batch
    from ddrug.models import Drug, MIC_COADD
    # from dsummary.utils.upd_sum_cmpbatch import sum_structure_sc
    # from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    # from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from adjcoadd.constants import COMPOUND_SEP
    # from django.db.models import Q

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'MIC_COADD':
        if prgArgs.runid:
            qryMIC = MIC_COADD.objects.filter(run_id=prgArgs.runid)
        elif prgArgs.new:
            qryMIC = MIC_COADD.objects.filter(bp_profile='')
        elif prgArgs.overwrite:
            qryMIC = MIC_COADD.objects.all()
        else:
            qryMIC = None
            
        nMIC = qryMIC.count()
        logger.info(f" [{prgArgs.table}] {nMIC} for {prgArgs.runid}")
        if qryMIC:

            OutNumbers = {'Processed':0,'New':0, 'Uploaded':0,'Empty':0, 'Failed':0}

            for mic in tqdm(qryMIC, total=nMIC, desc=prgArgs.table):
                validStatus = True 
                OutNumbers['Processed'] += 1
                mic.calc_breakpoint()

                mic.set_defaults_model()
                validDict = mic.validate_fields()
                if validDict:
                    validStatus = False
                    OutNumbers['Failed'] += 1 
                    logger.error(f"[{prgArgs.table}] {validDict}")

                if prgArgs.upload and validStatus:
                    mic.save()
                    OutNumbers['Uploaded'] += 1

        logger.info(f" [{prgArgs.table}] {OutNumbers}")





#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [CompoundID]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

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