import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from functools import reduce
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django
from oraABase.oraABase import openABase, get_CompoundBatch
from rdkit import Chem

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "RegChem_ABase"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def get_AbaseStructures(test=0):

    renameCol = {
        "labware_addon":  "labware_notes",
        "material":       "plate_material",
        "work_volume" :   "working_volume",
    }

    replaceValues = {
      'plate_size':{'384w':384,'96w':96,},
    }

    _SQL = "Select ObjdID, ObjsMolFormula, ObjsMolMassValue, ObjsMolFile  from ChemStruct "

    if test>0:
        _SQL += f" Fetch First {test} Rows Only "

    ABaseDB = openABase()
    logger.info(f"[ChemStructure] ... ")
    _DF = pd.DataFrame(ABaseDB.get_dict_list(_SQL))
    nTotal = len(_DF)
    logger.info(f"[ChemStructure] {nTotal} ")
    ABaseDB.close()

    # logger.info(f"DF - Rename Columns {len(renameCol)}")
    # _DF.rename(columns=renameCol, inplace=True)

    # logger.info(f"DF - Replace Values {len(replaceValues)}")
    # for k in replaceValues:
    #     _DF[k].replace(replaceValues[k],inplace=True)

    return(_DF)


def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dsample.models import Library, Library_Compound, Compound_Batch
    from dchem.models import Chem_Structure,Chem_Salt
    from dchem.utils.mol_std import get_Structure_Type, get_MF_Smiles, SaltDict_to_SaltCode, Smiles_to_Mol, SaltDictList_to_SaltCode
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # ABase ChemStructure -------------------------------------------------------------
    if prgArgs.table == "ChemStructure" :

        OutName = "[ChemStructure]"
        OutDict = []
        OutFile = f"UpdateABaseStruct_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}



        strSQL = "Select ObjdID, ObjsMolFormula, ObjsMolMassValue, ObjsMolFile  from ChemStruct "
        if int(prgArgs.test)>0:
            strSQL += f" Fetch First {int(prgArgs.test)} Rows Only "

        ABaseDB = openABase()

        print(f"{OutName} ---------------------------------------------------------")
        logger.info(f"[ChemStructure] ... ")
        strDF = pd.DataFrame(ABaseDB.get_dict_list(strSQL))
        nTotal = len(strDF)
        logger.info(f"[ChemStructure] {nTotal} ")
        print("--------------------------------------------------------------------")
        print(f"{OutName} {strDF.columns} ")

        for idx,row in tqdm(strDF.iterrows(),  total=len(strDF), desc="ABase Structure"):
            OutNumbers['Processed'] += 1
            if row['objsmolfile']:
                _molblock = row['objsmolfile'].read()
                smol = Chem.MolFromMolBlock(_molblock)
                if smol:
                    smol = Chem.MolFromMolBlock(_molblock,sanitize=False)
                    i = 1
                    #print(f" [{row['objdid']}] {Chem.Descriptors.MolWt(smol)} {row['objsmolmassvalue']}")
                if not smol:
                    print(f" [{row['objdid']}]  Unable to convert Structure ")
            # else:
            #     print(f" [{row['objdid']}]  No Structure ")

        ABaseDB.close()

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='regChem ABase', 
                                description="Register MCC structures from ABase")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [CompoundID]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.django == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.django == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.django == 'Laptop':
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================