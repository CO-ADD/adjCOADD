#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
from oraCastDB.oraCastDB import openCastDB
from pgCastDB.pgCastDB import openCoaddDB

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Check_AssayData"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------


#def convert_oraCmpBatch_djCmpBatch(CmpLst):

def main(prgArgs):

    if 'AssayData' in prgArgs.table :

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"checkAssayData_inORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        oraDB = openCastDB()
        djDB = openCoaddDB()

        AssayTable = {
            "AssayData_MIC": ['assaydata_mic','dscreen.assaydata_mic']
        }

        if prgArgs.table in AssayTable:
        # by oraCastDB
            oraSQL = f"Select testplate_id, testwell_id From {AssayTable[prgArgs.table][0]} where is_migrated < 1"
            nEntries = oraDB.nCount(f"Select count(1) From {AssayTable[prgArgs.table][0]} where is_migrated < 1" )
            logger.info(f"{OutName} {nEntries} ")

            oraDB.exec(oraSQL)  
            sql_columns = [i[0].lower() for i in oraDB.cursor.description]
            logger.info(sql_columns)

            updList = []
            for crow in tqdm(oraDB.cursor, total=nEntries, desc="[Check Assay]"):
                row = dict()
                for col in sql_columns:
                    row[col.lower()] = crow[sql_columns.index(col)]

                updVal = -1

                _ncount = djDB.nCount(f"Select count(1) From {AssayTable[prgArgs.table][1]} where testplate_id = '{row['testplate_id']}' and testwell_id = '{row['testwell_id']}' ")

                if _ncount > 0 :
                    updList.append((1,row['testplate_id'],row['testwell_id']))
                    
            logger.info(f"{OutName} {nEntries} -> {len(updList)}")
            if len(updList)>0 and prgArgs.upload:
                for row in tqdm(updList, desc="[Update oraProjects]"):
                    oraDB.exec(f"Update {AssayTable[prgArgs.table][1]} Set is_migrated = {row[0]} Where testplate_id = '{row[1]}' and testwell_id = '{row[2]}'",commit=True)     

        oraDB.close()
        djDB.close()

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

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    # prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    # prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    main(prgArgs)
    print("-------------------------------------------------------------------")

#==============================================================================