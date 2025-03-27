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

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_SMIcsv"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    from applib.data.set_fielddata import set_model_arrayfields, set_dictFields, set_model_dicts
    from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock, OrgBatch_Image

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # Table -------------------------------------------------------------
    if prgArgs.table == "Stock":


        ExcelFile = "C:\Data\LMIC\OrgDB_Curation.xlsx"
        ExcelSheet = "Stock_24Oct2024"

        stock_df = pd.read_excel(ExcelFile,sheet_name = ExcelSheet)

        rmColumns = ['OrganismID','BatchID','ID','X','Passages']

        appuser = ApplicationUser.get(prgArgs.appuser)
        empty_date = datetime.date(2009, 1, 1)

        for idx,row in tqdm(stock_df.iterrows(),total=len(stock_df)):

            if row['n_left'] > 0:

                #print(row)
                for rmCol in rmColumns:
                    if rmCol in row:
                        del row[rmCol]

                djBiologist = ApplicationUser.get(row['biologist'])
                djOrgBatch = Organism_Batch.get(row['orgbatch_id'])

                if djOrgBatch is not None:
                    djStock =  OrgBatch_Stock.get(None,OrgBatchID=djOrgBatch,StockDate=row['stock_date'],StockType=row['stock_type'])
                    if djStock is None:
                        djStock = OrgBatch_Stock()
                        djStock.orgbatch_id = djOrgBatch
                    row.pop('orgbatch_id')
 
                    djStock.stock_type = Dictionary.get(djStock.DICTIONARY_FIELDS["stock_type"],row['stock_type'])
                    row.pop('stock_type')

                    djStock.biologist = djBiologist
                    row.pop('biologist')

                    #print(djStock.stock_date)
                    #print(f" {djOrgBatch} {row['stock_date']} {row['stock_note']}")
                    if row['stock_date'] is pd.NaT:
                        row['stock_date'] = empty_date

                    # location_rack, location_column, location_slot => str(int())  

                    # set values in instance
                    for e in row.to_dict():
                        setattr(djStock,e,row[e])


                    djStock.set_defaults_model()
                    validDict = djStock.validate_fields()

                    if validDict:
                        logger.info(f" XX {djStock} {validDict} ")
                    else:
                        # --- Upload ---------------------------------------------------------
                        if prgArgs.upload:
                            djStock.save(user=appuser)
                else:
                    print(f" [OrgBatchID] {row['orgbatch_id']} not found {djBiologist}")

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
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
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
