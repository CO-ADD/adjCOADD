#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from oraCastDB.oraCastDB import openCastDB

from tqdm import tqdm
# from zUtils import zData

import django

# Logger ----------------------------------------------------------------
import logging
# logTime= datetime.datetime.now()
# logName = "Upload_VitekPDF"
# #logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
# logLevel = logging.INFO 

logger = logging.getLogger(__name__)
# logging.basicConfig(
#     format="[%(name)-20s] %(message)s ",
# #    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#     handlers=[logging.StreamHandler()],
#     level=logLevel)

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "Upload_VitekPDFDD"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
        #handlers=[logging.StreamHandler()],
        level=logLevel)
    #-----------------------------------------------------------------------------


    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from applib.data.str_lists import join_lst
    from dorganism.utils.utils import reformat_OrgBatchID

    from ddrug.models import Drug, MIC_COADD
    from dorganism.models import Organism_Batch
    from dscreen.models import Screen_Run

    import ddrug.utils.vitek as Vitek
    #from ddrug.utils.import_drug import imp_VitekCard_fromDict
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "Vitek" :
        print("--> imp Vitek  ---------------------------------------------------------")

        dirName = ""
        fileNames = []
        if prgArgs.file:
            if os.path.isfile(prgArgs.file):
                dirName,fileName = os.path.split(prgArgs.file)
                fileNames = [fileName]
        elif prgArgs.directory:
            if os.path.isdir(prgArgs.directory):
                dirName = prgArgs.directory
                files = os.listdir(dirName)
                fileNames = [f for f in files if f.endswith(".pdf")]

        logger.info(f"[Upd_djCOADD] {fileNames} from folder {dirName} [Upload: {prgArgs.upload}]") 
        if len(fileNames) > 0 :
            for fileName in fileNames:
                Vitek.upload_VitekPDF(dirName,fileName,OrgBatchID=None,upload=prgArgs.upload) 
                #lCards,lID,lAST = Vitek.process_VitekPDF(dirName,fileName)



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
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of rows to test")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
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
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
