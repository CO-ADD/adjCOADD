#
#
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
    prgParser = argparse.ArgumentParser(prog='change_OrgDB_Data', 
                                description="Changing data to adjCOADD from Excel")
    # prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    # prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    # prgParser.add_argument("-o","--orgbatch",default=None,required=False, dest="orgbatch", action='store', help="OrganismBatch ID")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    # prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD"
    uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"

    xlFiles = {
        'Application': "ApplicationData_v05.xlsx",
        'Drug': "DrugData_v04.xlsx",
        'MIC': "LMIC_Data_v06.xlsx",
        'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser
    import b_change_dOrganism as dOrg

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "ChangeOrgDB"
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

    appuser = ApplicationUser.get(prgArgs.appuser)

    #if prgArgs.file:
    #    dOrg.rename_OrgName_xls(prgArgs.file,XlsSheet="New OrgName",
    #                            upload=prgArgs.upload,uploaduser=appuser)

    dOrg.rename_OrgID('GN_0696','Bacillus velezensis',None,{'Old Strain Code':'Kp BRf 79136','strain_code':'Ba.ve BRf 79136'})

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
