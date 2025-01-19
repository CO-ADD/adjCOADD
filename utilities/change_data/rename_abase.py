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
logName = "Update_ABase_MCC"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)



def rename_cmpbatchlst(cmpbatch_lst):
    isChanged = False
    nCmpBatchLst = []
    for cmpbatch in cmpbatch_lst:
        if 'MCC0' in cmpbatch:
            nCmpBatchLst.append(cmpbatch.replace('MCC0','MCC_0'))
            isChanged = True
        else:
            nCmpBatchLst.append(cmpbatch)
    return(nCmpBatchLst,isChanged)

def move_cmpbatch_instance(djInst,upload=False):
    nCmpBatchLst,isChanged = rename_cmpbatchlst(djInst.cmpbatch_lst)
    nChg = 0
    nUpl = 0
    if isChanged:
        nChg = 1
        if upload:
            nUpl = 1
            # djInst.astatus = -9
            # djInst.save()
            djInst.cmpbatch_lst = nCmpBatchLst
            djInst.astatus = 0
            djInst.save()
    return(nChg,nUpl)

#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    from dscreen.models import AssayData_CC50,AssayData_HC50,AssayData_MIC

    OutNumbers = {'Processed':0, 'Empty':0, 'Failed':0, 'New':0, 'Uploaded':0,}
    nEntries = 0
    # ----------------------------------------------------------------------------------
    if prgArgs.table in 'AssayData_HC50':
        qry = AssayData_HC50.objects.all()
        nEntries = qry.count()
    if prgArgs.table in 'AssayData_MIC':
        qry = AssayData_MIC.objects.all()
        nEntries = qry.count()
        qry = AssayData_MIC.objects.all().iterator(chunk_size=100)

    logger.info(f" [{prgArgs.table}] : {nEntries}")
    # ----------------------------------------------------------------------------------
    if nEntries > 0:     
        for djInst in tqdm(qry, total=nEntries):
            OutNumbers['Processed'] += 1
            nChg,nUpl = move_cmpbatch_instance(djInst,upload=prgArgs.upload)
            OutNumbers['Uploaded'] += nUpl
            OutNumbers['New'] += nChg
        logger.info(f" [{prgArgs.table}] : {OutNumbers}")
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
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.django == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.django == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.django == 'Laptop':
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================    