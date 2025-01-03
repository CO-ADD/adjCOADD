##
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
from oraCastDB.oraCastDB import openCastDB
from pgCastDB.pgCastDB import openCoaddDB

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_Direct"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)


def main(prgArgs):
    oraDB = openCastDB()
    djDB = openCoaddDB()
    if prgArgs.table == "TestWells" :

        OutName = "[TestWells]"
        oraSQL = "Select plate_id, well_id, active From TestWell "

        if int(prgArgs.test) > 0:
            oraSQL += f" Fetch First {int(prgArgs.test)} Rows Only "
            nWells = int(prgArgs.test)
        else:
            nWells = oraDB.nCount(f"Select count(1) From TestWell " )
        logger.info(f"{OutName} {nWells} ")


        oraDB.exec(oraSQL)  
        sql_columns = [i[0].lower() for i in oraDB.cursor.description]
        logger.info(sql_columns)

        for crow in tqdm(oraDB.cursor, total=nWells, desc=OutName):

            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]

            djSQL = f"Update dplate.testwell Set act_type = '{row['active']}' Where plate_id = '{row['plate_id']}' and well_id = '{row['well_id']}' "
            if prgArgs.upload:
                djDB.exec(djSQL,commit=True) 

    djDB.close()
    oraDB.close()


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
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    # prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    # prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    main(prgArgs)
    print("-------------------------------------------------------------------")
#==============================================================================
