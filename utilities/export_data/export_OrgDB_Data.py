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
    djDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD"
    outputDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/export_data/Output"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        outputDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/export_data/Output"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        outputDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/export_data/Output"

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    import d_export_Vitek as dVitek

   # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "ExportOrgDB"
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
    OutBase = "Output"
    #uploadFile = os.path.join(uploadDir,"ApplicationData_2022_11_21_JZuegg_v01.xlsx")

    choiceTables = ['OrgBatch','Vitek','MIC',
                    'wgsFastA','wgsAMR',
                    ]
    if prgArgs.table in choiceTables:
        logger.info(f"[Exp_djCOADD] Table: {prgArgs.table}") 

        if prgArgs.table == 'OrgBatch':
            OutDir = os.path.join(OutBase,"OrgBatch")
            dVitek.export_OrgBatch(OutDir)

        if prgArgs.table == 'Vitek':
            OutDir = os.path.join(OutBase,"Vitek")
            dVitek.export_Vitek(OutDir)

        if prgArgs.table == 'MIC':
            OutDir = os.path.join(OutBase,"Antibiogram")
            dVitek.export_Antibiogram(OutDir)

        if prgArgs.table == 'wgsFastA':
            OutDir = os.path.join(OutBase,"FastA")
            dVitek.export_wgsFastA(OutDir)

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
