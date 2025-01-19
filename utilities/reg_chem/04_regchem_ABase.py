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
    from dsample.models import ABase_Compound
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # ABase ChemStructure -------------------------------------------------------------
    if prgArgs.table == "Update_ChemStructure" :

        OutName = "[Update_ChemStructure]"
        OutDict = []
        OutFile = f"UpdateABaseStruct_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New':0, 'Uploaded':0, 'Failed':0}

        strSQL = "Select ObjdID, ObjsMolFormula, ObjsMolMassValue, ObjsMolFile  from ChemStruct "
        if int(prgArgs.test)>0:
            strSQL += f" Fetch First {int(prgArgs.test)} Rows Only "

        ABaseDB = openABase()

        print(f"{OutName} ---------------------------------------------------------")
        logger.info(f"[ChemStructureReg] ... ")
        strDF = pd.DataFrame(ABaseDB.get_dict_list(strSQL))
        nTotal = len(strDF)
        logger.info(f"[ChemStructureReg] {nTotal} ")
        print("--------------------------------------------------------------------")
        print(f"{OutName} {strDF.columns} ")

        for idx,row in tqdm(strDF.iterrows(),  total=len(strDF), desc="ABase Structure"):
            NewEntry = False
            validStatus = True
            OutNumbers['Processed'] += 1
            djCmp =  ABase_Compound.get(row['objdid'])
            if djCmp is None:
                NewEntry = True
                djCmp = ABase_Compound()
                djCmp.compound_id = row['objdid']
                OutNumbers['New'] += 1
            
            if row['objsmolfile']:
                _molblock = row['objsmolfile'].read()
                djCmp.reg_molfile = _molblock
                djCmp.reg_mw = row['objsmolmassvalue']
                djCmp.reg_mf = row['objsmolformula']

            djCmp.init_fields()
            validDict = djCmp.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    logger.warning('Warning',k,validDict[k],'-')
                OutDict.append(row)

            if validStatus:
                if prgArgs.upload:
                    if NewEntry or prgArgs.overwrite:
                        OutNumbers['Uploaded'] += 1
                        djCmp.save(user=prgArgs.appuser)
        ABaseDB.close()
    
    elif prgArgs.table == "Reg_ChemStructure":


        OutName = "[Reg_ChemStructure]"
        OutDict = []
        OutFile = f"RegABaseStruct_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New':0, 'Uploaded':0, 'Failed':0,'Empty':0, 'not STD':0}


        print(f" {OutName} ---------------------------------------------------------")
        logger.info(f"[Reg_ChemStructure] ... ")
        if int(prgArgs.test)>0:
            qryStr = ABase_Compound.objects.all()[:int(prgArgs.test)]
            nEntry = int(prgArgs.test)
        else:
            qryStr = ABase_Compound.objects.all()
            nEntry = qryStr.count()
        logger.info(f" {OutName}: {nEntry} ")

        for djMCC in tqdm(qryStr, total=nEntry, desc='[CmpBatcheLsts]'):
            OutNumbers['Processed'] += 1

            if djMCC.reg_molfile:
                _mol = Chem.MolFromMolBlock(djMCC.reg_molfile)
                if _mol:
                    validStatus = True
                    Chem.Kekulize(_mol)
                    _MolType,_Metal,_IsMet = get_Structure_Type(_mol,None)
                    djMCC.structure_type = _MolType
                    djMCC.structure_metal = _Metal

                    _smi = Chem.MolToSmiles(_mol, canonical=True, isomericSmiles=True, kekuleSmiles=True)

                    djChem = Chem_Structure.get_exact(_mol)
                    if djChem is None:
                        djChem = Chem_Structure()
                        djChem.smol = _mol
                        OutNumbers['New'] += 1

                    djChem.init_fields()
                    validDict = djChem.validate_fields()
                    
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            logger.warning(f"{k}: {validDict[k]}")
                            
                    if prgArgs.upload and validStatus:
                        djChem.save()
                        _csid = djChem.structure_id
                        OutNumbers['Uploaded'] += 1 
                        # Reload to get MW
                        djChem.get(_csid)
                        djMCC.structure_id = djChem

                    
                else:
                    _mol = Chem.MolFromMolBlock(djMCC.reg_molfile, sanitize=False)
                    if _mol:
                        OutNumbers['not STD'] += 1
                        _MolType,_Metal,_IsMet = get_Structure_Type(_mol,None)
                        djMCC.structure_type = _MolType
                        djMCC.structure_metal = _Metal
                    else:    
                        OutNumbers['Failed'] += 1

                validStatus = True
                validDict = djMCC.validate_fields()
                djMCC.init_fields()
                validDict = djMCC.validate_fields()
                
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")

                if prgArgs.upload and validStatus:
                    djMCC.save()

            else:
                OutNumbers['Empty'] += 1
        logger.info(f"[{prgArgs.table}] {OutNumbers}")

# In [14]: m = Chem.MolFromSmiles('F[P-](F)(F)(F)(F)F.CN
# (C)C(F)=[N+](C)C',sanitize=False)

# In [15]: m.UpdatePropertyCache(strict=False)

# In [16]:
# Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
# Out[16]: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE



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