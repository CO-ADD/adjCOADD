#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from zChem.zMolStandardize import SmiStandardizer_DB
from rdkit import Chem

from tqdm import tqdm
# from zUtils import zData

import django
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    # Logger ----------------------------------------------------------------
    import logging
    logTime= datetime.datetime.now()
    logName = "regChem_01cUpdChem_COADD"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
        #handlers=[logging.StreamHandler()],
        level=logLevel)
    #-----------------------------------------------------------------------------


    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from apputil.utils.data import Dict_to_StrList
    from dsample.models import Project, COADD_Compound, Sample, Convert_ProjectID, Convert_CompoundID
    from dchem.models import Chem_Salt
    from dchem.utils.mol_std import get_Structure_Type_Smiles, get_MF_Smiles, SaltDict_to_SaltCode

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "COADD_Compound" :


        MolStd = SmiStandardizer_DB(chemdb=Chem_Salt) 

        logger.info("--> COADD_Compound ---------------------------------------------------------")
        #qryCmpd = COADD_Compound.objects.exclude(reg_smiles="")
        logger.info(f"[CO-ADD Compound] Reading {prgArgs.file}")
        dfCmpd = pd.read_excel(prgArgs.file)
        nCmpd = len(dfCmpd)   
        logger.info(f"[CO-ADD Compound] {nCmpd}")
        logger.info("-------------------------------------------------------------------------")
        OutFile = f"updChem_COADD_{logTime:%Y%m%d_%H%M%S}.xlsx"

        outNumbers = {'Proc':0,'Updated Compounds':0}

        for idx,row in tqdm(dfCmpd.iterrows(), total=nCmpd, desc="Updating Compounds"):
            
            updFields = ['std_status','std_process','std_smiles','std_salt','std_mw_extra','std_smiles_extra','std_issues']
            updDict = ['compound_type']
            djCmpd = COADD_Compound.get(row['compound_id'])
            if djCmpd:
                validStatus = True
                set_dictFields(djCmpd,row,updFields)
                # set_arrayFields(djCmpd,row,arrayFields)
                set_Dictionaries(djCmpd,row,updDict)
                
                djCmpd.clean_Fields()
                validDict = djCmpd.validate()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")

                if validStatus and prgArgs.upload:
                    outNumbers['Updated Compounds'] += 1
                    djCmpd.save()

        logger.info(f"[CO-ADD Compound] {outNumbers}")


    
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
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of rows to test")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=True, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

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
