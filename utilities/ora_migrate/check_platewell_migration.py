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
logName = "Check_WellMigration"
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

    #-----------------------------------------------------------------------------
    if prgArgs.table == "TestWells" :

        OutName = "[TestWells]"
        OutDict = []
        OutFile = f"checkTestWells_inORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        oraDB = openCastDB()
        djDB = openCoaddDB()

        djSQL = "Select plate_id, well_id from dplate.testwell "
        nWells = djDB.nCount("Select count(1) From dplate.testwell" )
        logger.info(f"{OutName} {nWells} ")


        djDB.exec(djSQL)  
        sql_columns = [i[0].lower() for i in djDB.cursor.description]
        logger.info(sql_columns)


        
        for crow in tqdm(djDB.cursor, total=nWells, desc=OutName):
            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]
            updSQL = f"Update TestWell Set is_migrated = 1 Where Plate_ID = '{row['plate_id']}' and Well_ID = '{row['well_id']}' "
            oraDB.exec(updSQL,commit=True)

        oraDB.close()
        djDB.close()

    #-----------------------------------------------------------------------------
    if prgArgs.table == "MasterWells" :

        OutName = "[MasterWells]"
        OutDict = []
        OutFile = f"checkMasterWells_inORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        oraDB = openCastDB()
        djDB = openCoaddDB()

        djSQL = "Select plate_id, well_id, barcode from dplate.masterwell"
        nWells = djDB.nCount("Select count(1) From dplate.masterwell" )
        logger.info(f"{OutName} {nWells} ")


        djDB.exec(djSQL)  
        sql_columns = [i[0].lower() for i in djDB.cursor.description]
        logger.info(sql_columns)
        
        for crow in tqdm(djDB.cursor, total=nWells, desc=OutName):
            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]
            updSQL = f"Update MasterWell Set is_migrated = 1 Where Plate_ID = '{row['plate_id']}' and Well_ID = '{row['well_id']}' "
            oraDB.exec(updSQL,commit=True)

    #-----------------------------------------------------------------------------
    if prgArgs.table == "Barcodes" :

        OutName = "[Barcodes]"
        OutDict = []
        OutFile = f"checkBarcodes_inORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'Found':0, 'Missing':0, 'Duplicates':0}

        oraDB = openCastDB()
        djDB = openCoaddDB()
        
        oraSQL = "Select plate_id, well_id, barcode, compound_id, conc, conc_unit From MasterWell Where barcode is not Null"

        # djSQL = "Select plate_id, well_id, barcode from dplate.masterwell"
        nBarcodes = oraDB.nCount("Select count(1) from MasterWell Where barcode is not Null" )
        logger.info(f"{OutName} {nBarcodes} in oraCastDB")

        oraDB.exec(oraSQL)  
        sql_columns = [i[0].lower() for i in oraDB.cursor.description]
        logger.info(sql_columns)
        
        for crow in tqdm(oraDB.cursor, total=nBarcodes, desc=OutName):
            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]
            djSQL = f"Select count(1) From  dplate.masterwell Where barcode = '{row['barcode']}' "
            n_barcodes = djDB.nCount(djSQL)
            
            if n_barcodes == 0:
                OutNumbers['Missing'] += 1 
                logger.warning(f" Barcode not found in djCastDB {row['plate_id']} {row['well_id']} {row['barcode']} {row['compound_id']} {row['conc']} {row['conc_unit']}")
            elif n_barcodes == 1:
                OutNumbers['Found'] += 1 
            elif n_barcodes > 1:
                OutNumbers['Duplicates'] += 1 
        oraDB.close()
        djDB.close()

        logger.info(f" [{OutName}] {OutNumbers}")

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [Barcodes/MasterWells/TestWells]")
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

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    main(prgArgs)
    print("-------------------------------------------------------------------")

#==============================================================================