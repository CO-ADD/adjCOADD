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
logName = "Check_Compounds"
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

    if prgArgs.table == "Compounds" :

        OutName = "[Compounds]"
        OutDict = []
        OutFile = f"checkCompounds_inORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        oraDB = openCastDB()
        djDB = openCoaddDB()

        # by oraCastDB
        oraSQL = "Select compound_id From Compound where is_migrated < 1"
        nWells = oraDB.nCount("Select count(1) From Compound where is_migrated < 1" )
        logger.info(f"{OutName} {nWells} ")

        oraDB.exec(oraSQL)  
        sql_columns = [i[0].lower() for i in oraDB.cursor.description]
        logger.info(sql_columns)

        updList = []
        for crow in tqdm(oraDB.cursor, total=nWells, desc="[Check Compounds]"):
            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]
            oraCID = crow[0]

            updVal = -1

            djDB.exec(f"Select compound_id From dsample.convert_compoundid where ora_compound_id = '{oraCID}' ")
            entry = djDB.cursor.fetchone()
            if entry:
                convID = entry[0]
            else:
                convID = None

            if convID:
                djDB.exec(f"Select compound_id From dsample.coadd_compound where compound_id = '{convID}' ")
                entry = djDB.cursor.fetchone()
                if entry:
                    coaddID = entry[0]
                else:
                    coaddID = None

                djDB.exec(f"Select cmpbatch_id From dsample.cmpbatch where cmpbatch_id = '{convID}' ")
                entry = djDB.cursor.fetchone()
                if entry:
                    cmpbatchID = entry[0]
                else:
                    cmpbatchID = None

                if cmpbatchID:
                    updVal = 0
                    if coaddID:
                        updVal = 1

            if updVal > -1:
                updList.append((oraCID,updVal))

        if len(updList)>0:
            for row in tqdm(updList, desc="[Update oraCompounds]"):
                oraDB.exec(f"Update Compound Set is_migrated = {row[1]} Where Compound_ID = '{row[0]}' ",commit=True)     

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
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [Compounds]")
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