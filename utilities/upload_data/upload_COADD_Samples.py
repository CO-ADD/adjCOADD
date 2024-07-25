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
tqdm.pandas()
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
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "D:/Code/zdjCode/adjCOADD"
    # uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    # if prgArgs.database == 'Work':
    #     djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # elif prgArgs.database == 'WorkLinux':
    #     djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"

    # xlFiles = {
    #     'Application': "ApplicationData_v05.xlsx",
    #     'Drug': "DrugData_v04.xlsx",
    #     'MIC': "LMIC_Data_v06.xlsx",
    #     'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    # }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dsample.models import Project, COADD_Sample    
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "UploadOrgDB"
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

    choiceTables = ['COADD_Sample',
                    ]
    if prgArgs.table in choiceTables:

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

        if prgArgs.table == 'COADD_Sample' and prgArgs.excel:
            print(f"Reading {prgArgs.excel}.[SampleID]")
            df = pd.read_excel(prgArgs.excel,sheet_name='SampleID')
            df.compound_name = df.compound_name.fillna('')
            print(df.columns)
            outDict = []    
            for idx,row in tqdm(df.iterrows(), total=df.shape[0]):
                if row['stype'] in ['C0','CX']:
                    #print(f"{row['old_compound_id']} {row['sample_id']}")                
                    djSmp = COADD_Sample.get(row['sample_id'])
                    if djSmp is None: 
                        djSmp = COADD_Sample()
                        djPrj = Project.get(row['project_id'])
                        if djPrj is not None: 
                            djSmp.project_id = djPrj
                            djSmp.old_compound_id = row['old_compound_id']
                            djSmp.sample_id = row['sample_id']
                            djSmp.sample_code = row['compound_code']
                            djSmp.sample_name = row['compound_name']
                            if prgArgs.upload:
                                djSmp.clean_Fields()
                                djSmp.save()
                        else:
                            row['Issue'] = f"No Project ID"
                            #row['IssueValue'] = f"{row['project_id']}"
                            outDict.append(row)
                    else:
                        row['Issue'] = f"Sample_ID Exists"
                        #row['IssueValue'] = f"{row['sample_id']}; {row['old_compound_id']}"
                        outDict.append(row)
                else:
                    row['Issue'] = f"Not C0 or CX"
                    #row['IssueValue'] = f"{row['sample_id']}; {row['old_compound_id']}"
                    outDict.append(row)
                                                    
            outDF = pd.DataFrame(outDict)
            outDF.to_excel('UploadSamples_Issues.xlsx')

    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
