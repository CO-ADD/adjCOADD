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

import logging
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "D:/Code/zdjCode/adjCOADD"

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dsample.models import Project
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "Upload_ConvertID"
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

   # Table -------------------------------------------------------------
    runTables = ["CompoundID","ProjectID"]

    if "ProjectID" in runTables:
        ExcelFile = "D:/Upload/CastDB/Old2New_ID.xlsx"
        SheetName = "ProjectID"
        print(f"[Reading Excel] {ExcelFile} {SheetName} ")
        prjDF = pd.read_excel(ExcelFile,sheet_name=SheetName)
        prjDF.project_name = prjDF.project_name.fillna('')
        
        for idx,row in tqdm(prjDF.iterrows(), total=prjDF.shape[0]):
            djE = Convert_ProjectID.get(row['ora_project_id'])
            if djE is None:
                djE = Convert_ProjectID()
                djE.ora_project_id = row['ora_project_id']
                djE.project_id = row['project_id']
                djE.project_name = row['project_name']
                if prgArgs.upload:
                    djE.clean_Fields()
                    djE.save()
            else:
                print(f"[Exists already] {row['ora_project_id']} {row['project_id']} ")

    if "CompoundID" in runTables:
        ExcelFile = "D:/Upload/CastDB/Old2New_ID.xlsx"
        SheetName = "CompoundID"
        print(f"[Reading Excel] {ExcelFile} {SheetName} ")
        cmpDF = pd.read_excel(ExcelFile,sheet_name=SheetName)
        cmpDF.compound_name = cmpDF.compound_name.fillna('')
        cmpDF.compound_code = cmpDF.compound_code.fillna('')
        
        for idx,row in tqdm(cmpDF.iterrows(), total=cmpDF.shape[0]):
            djE = Convert_CompoundID.get(row['ora_compound_id'])
            if djE is None:
                djE = Convert_CompoundID()
                djE.ora_compound_id = row['ora_compound_id']
                djE.compound_id = row['compound_id']
                djE.compound_name = row['compound_name']
                djE.compound_code = row['compound_code']
                djE.sample_type = row['stype']
                djE.project_id = row['project_id']
                if prgArgs.upload:
                    djE.clean_Fields()
                    djE.save()
            else:
                print(f"[Exists already] {row['ora_compound_id']} {row['compound_id']} ")
    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================