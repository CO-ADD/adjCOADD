#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from zSql import zSqlConnector
from rdkit import Chem 

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


def openChemDB(User='chemdb', Passwd='CHEMDB'):
    db = zSqlConnector.Oracle()
    logger.info("[ChemDB] Connecting to oraChem ")
    db.open(User,Passwd,"imb-coadd-db.imb.uq.edu.au","1521","coadb")
    return(db)


def get_oarSalt(test=0):

    saltSQL = """
            Select SALT_ID, SALT_NAME, SALT_TYPE, SMILES,
                INCHIKEY,
                MW,MF,MASS,NATOM,CHARGE,H_EQUIV
            From  ChemSalt
            """
    
    if test>0:
        saltSQL += f" Fetch First {test} Rows Only "

    oraChemDB = openChemDB()
    logger.info(f"[ChemDB] ... ")
    cmpdLst = oraChemDB.get_dict_list(saltSQL)
    nTotal = len(cmpdLst)
    logger.info(f"[ChemDB]  {nTotal} Salts")
    oraChemDB.close()

    return(cmpdLst)



#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from dchem.models import Chem_Structure, Chem_Salt
    from dsample.models import Library, Library_Compound

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # Table -------------------------------------------------------------
    if prgArgs.table == "Salt":

        cpyFields = ['compound_name',
                    'reg_smiles',
                    ]
        
        appuser = ApplicationUser.get(prgArgs.appuser)


        cmpdLst = get_oarSalt(int(prgArgs.test))


        outNumbers = {'Proc':0,'New Compounds':0,'Upload Compounds':0, 'New Samples': 0, 'Upload Samples': 0, 'Failed': 0}
        outDict = []
        new_compound = False
        
        for row in tqdm(cmpdLst):
            outNumbers['Proc'] += 1
            djSalt = Chem_Salt.get(row['salt_id'])
            if not djSalt:
                djSalt = Chem_Salt()
                djSalt.salt_id = row['salt_id']
                new_compound = True

            set_dictFields(djSalt,row,['salt_name','smiles','mw','mf','h_equiv'])
            # set_arrayFields(djCmpd,row,arrayFields)
            set_Dictionaries(djSalt,row,['salt_type'])

            validStatus = True

            djSalt.clean_Fields()
            validDict = djSalt.validate()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('Warning',k,validDict[k],'-')
                outDict.append(row)

            if validStatus:
                if prgArgs.upload:
                    if new_compound or prgArgs.overwrite:
                        outNumbers['Upload Compounds'] += 1
                        djSalt.save()

        print(f"[LibCompounds] :{outNumbers}")

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
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of rows to test")
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
