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


from zSql import zSqlConnector
from rdkit import Chem 

from tqdm import tqdm
# from zUtils import zData

from zDjango.djUtils import init_django_dir
import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ChEMBL"
logDir = "log"
logFileName = os.path.join(logDir,f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

if not os.path.isdir(logDir):
    os.mkdir(logDir)

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def openChEMBL(User='chembl', Passwd='chembl',DataBase="chembl",verbose=1):
    dbPG = zSqlConnector.PostgreSQL()
    dbPG.open(User,Passwd,"imb-coadd-db.imb.uq.edu.au",DataBase,verbose=verbose)
    return(dbPG)

#-----------------------------------------------------------------------------
def get_pgCompound(test=0):

    cmpdSQL = """
    Select s.canonical_smiles reg_smiles, s.molregno, m.chembl_id compound_code, m.pref_name compound_name, molecule_type compound_type
    From compound_structures s
     Left join molecule_dictionary m on s.molregno = m.molregno
    Where s.canonical_smiles is not Null
    """
    if test>0:
        cmpdSQL += f" Fetch First {test} Rows Only "

    pgCHEMBL = openChEMBL()
    logger.info(f"[ChEMBL] ... ")
    cmpdLst = pgCHEMBL.get_dict_list(cmpdSQL)
    nTotal = len(cmpdLst)
    logger.info(f"[ChEMBL]  {nTotal} Compounds")
    pgCHEMBL.close()

    return(cmpdLst)



#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir['djPrj'])
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    logging.getLogger().addHandler(logging.FileHandler(logFileName,mode='w'))

    from apputil.models import ApplicationUser, Dictionary
    from applib.data.set_fielddata import set_arrayFields, set_dictFields, set_Dictionaries
    from dchem.models import Chem_Structure, Chem_Salt
    from dsample.models import Library, Library_Compound

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    LibraryID = 'CHEMBL'

    # Table -------------------------------------------------------------
    if prgArgs.table == "Library":

        cpyFields = ['compound_name',
                    'reg_smiles',
                    ]
        
        appuser = ApplicationUser.get(prgArgs.appuser)
        djLib = Library.get(LibraryID)

        if djLib:

            cmpdLst = get_pgCompound(int(prgArgs.test))

            outNumbers = {'Proc':0,'New Compounds':0,'Upload Compounds':0, 'New Samples': 0, 'Upload Samples': 0, 'Failed': 0}
            outDict = []
            new_compound = False
            
            for row in tqdm(cmpdLst):
                djCmpd = Library_Compound.get(None,row['compound_code'],LibraryID)
                if not djCmpd:
                    djCmpd = Library_Compound()
                    djCmpd.compound_code = row['compound_code']
                    djCmpd.library_id = djLib
                    new_compound = True
                djCmpd.compound_desc = f"MolRegNo: {row['molregno']}"

                set_dictFields(djCmpd,row,['reg_smiles','compound_name'])
                # set_arrayFields(djCmpd,row,arrayFields)
                set_Dictionaries(djCmpd,row,['compound_type'])

                validStatus = True

                djCmpd.init_fields()
                validDict = djCmpd.validate_fields()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        print('Warning',k,validDict[k],'-')
                    outDict.append(row)

                if validStatus:
                    if prgArgs.upload:
                        if new_compound or prgArgs.overwrite:
                            outNumbers['Upload Compounds'] += 1
                            djCmpd.save()

            print(f"[{LibraryID}] :{outNumbers}")

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")

    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of rows to test")

    prgParser.add_argument("--django",default='Local',required=True, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
