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

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ConvertID"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from apputil.utils.data import Dict_to_StrList
    from dsample.models import Project, COADD_Compound, Sample, Convert_ProjectID, Convert_CompoundID
    from dchem.models import Chem_Salt
    from dchem.utils.mol_std import check_structure_type_smi, SaltDict_to_SaltCode

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "COADD_Compound" :


        MolStd = SmiStandardizer_DB(chemdb=Chem_Salt) 

        print("--> COADD_Compound ---------------------------------------------------------")
        #qryCmpd = COADD_Compound.objects.exclude(reg_smiles="")
        qryCmpd = COADD_Compound.objects.all()
        print(f"[CO-ADD Compound] {len(qryCmpd)}")
        print("-------------------------------------------------------------------------")
        OutFile = f"regChem_COADD_{logTime:%Y%m%d_%H%M%S}.xlsx"

        outNumbers = {'Proc':0,'Updated Compounds':0, 'Metal Compounds':0, 'Already Done': 0}

        for djCmpd in tqdm(qryCmpd):
            outNumbers['Proc'] += 1
            updated_sample = False

            # Check if this Standardisation has been done already 
            #if not djCmpd.std_status or djCmpd.std_status == 'Invalid' or prgArgs.overwrite:
            if djCmpd.reg_smiles or djCmpd.reg_mf:

                _MolType,_Metal,_IsMet = check_structure_type_smi(djCmpd.reg_smiles,djCmpd.reg_mf)
                djCmpd.std_structure_type = _MolType
                djCmpd.std_metal = _Metal
                updated_sample = True
                validStatus = True

            #if not djCmpd.std_status or djCmpd.std_status != 'Valid' or prgArgs.overwrite:

                # Excluded from SmiStandardizer - as molvs breaks any metal bonds
                # metal specific Standardizer is required, including OpenSmiles syntax for Metalcomplex
                if _IsMet==1:
                    djCmpd.std_status = 'Metal'
                    outNumbers['Metal Compounds'] += 1

                    validStatus = True
                    updated_sample = True

                # Non Metal complex structures
                elif djCmpd.reg_smiles and djCmpd.std_status != 'Valid':
 
                    _moldict, _saltdict, _iondict, _solvdict = MolStd.run_single(djCmpd.reg_smiles)

                    if _moldict['valid'] > 0:
                        djCmpd.std_status = 'Valid'
                        djCmpd.std_process = "Std"

                        djCmpd.std_smiles = _moldict['smi']
                        djCmpd.std_mw = _moldict['mw']
                        djCmpd.std_nfrag = _moldict['nfrag']

                        djCmpd.std_salt = SaltDict_to_SaltCode(_saltdict)
                        djCmpd.std_ion = SaltDict_to_SaltCode(_iondict)
                        djCmpd.std_solvent = SaltDict_to_SaltCode(_solvdict)
                        djCmpd.std_smiles_extra = _moldict['smiles_extra']
                        djCmpd.std_mw_extra = _moldict['mw_extra']

                        validStatus = True
                        updated_sample = True
                    else:
                        djCmpd.std_status = 'Invalid'
                        djCmpd.std_process = "Std"

                        validStatus = True
                        updated_sample = True

                djCmpd.clean_Fields()
                validDict = djCmpd.validate()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        print('Warning',k,validDict[k],'-')

                if validStatus and updated_sample and prgArgs.upload:
                    outNumbers['Updated Compounds'] += 1
                    djCmpd.save()
            else:
                outNumbers['Already Done'] += 1
        print(f"[CO-ADD Compound] {outNumbers}")


    
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
