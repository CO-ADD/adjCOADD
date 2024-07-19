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
    logName = "regChem_01eRegChem_COADD"
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
    from dchem.models import Chem_Structure, Chem_Salt
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
        qryCmpd = COADD_Compound.objects.filter(std_status='Valid')            
        nCmpd = qryCmpd.count()    
        logger.info(f"[CO-ADD Compound] {nCmpd} [Valid]")
        logger.info("-------------------------------------------------------------------------")
        OutFile = f"regChem_COADD_{logTime:%Y%m%d_%H%M%S}.xlsx"

        outNumbers = {'Proc':0,'Updated Compounds':0,
                      'New Samples':0, 'Updated Samples': 0,  
                      'New ChemStructures':0, 'Updated ChemStructures': 0}

        for djCmpd in tqdm(qryCmpd.iterator(), total=nCmpd, desc="Processing Compounds"):
            outNumbers['Proc'] += 1
            updated_sample = False

            # Check if this Standardisation has been done already 
            #if not djCmpd.std_status or djCmpd.std_status == 'Invalid' or prgArgs.overwrite:
            
            if djCmpd.std_nfrag == 1:
                updated_sample = True
                validStatus = True
                outNumbers['Updated Compounds'] += 1

                #------------------------------------------------------------
                djChem = Chem_Structure.get_bySmiles(djCmpd.std_smiles)
                if djChem is None:
                    djChem = Chem_Structure()
                    djChem.set_molecule(djCmpd.std_smiles)
                    djChem.nfrag = djCmpd.std_nfrag
                    outNumbers['New ChemStructures'] += 1
                    #logger.info(f"[CO-ADD Compound] New Chem_Structure {djCmpd.std_smiles}")


                    djChem.clean_Fields()
                    validDict = djChem.validate()
                    
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            logger.warning(f"{k}: {validDict[k]}")
                            
                    if prgArgs.upload and validStatus:
                        #djCmpd.std_process += ";ChemStructure"
                        djChem.save()
                        outNumbers['Updated ChemStructures'] += 1
                else: 
                    #logger.info(f"[CO-ADD Compound] Existing Chem_Structure {djChem}")
                #------------------------------------------------------------
                djSample = Sample.get(djCmpd.compound_id)
                if djSample is None:
                    djSample = Sample()
                    djSample.sample_id = djCmpd.compound_id
                    djSample.sample_source = 'COADD'
                    new_sample = True
                    outNumbers['New Samples'] += 1

                djSample.sample_code = djCmpd.compound_code
                djSample.structure_id = djChem
                djSample.structure_type = djCmpd.std_structure_type
                _salt_code = []
                if djCmpd.std_salt:
                    _salt_code.append(djCmpd.std_salt)   
                if djCmpd.std_ion:
                    _salt_code.append(djCmpd.std_ion)   
                if djCmpd.std_solvent:
                    _salt_code.append(djCmpd.std_solvent)                            
                djSample.salt_code = ";".join(_salt_code)
                
                djSample.smiles_extra = djCmpd.std_smiles_extra
                djSample.mw_extra = djCmpd.std_mw_extra
                djSample.full_mw = float(djSample.mw_extra) + float(djChem.mw)
                djSample.full_mf = get_MF_Smiles(djCmpd.std_smiles+djSample.smiles_extra)
                
                djSample.clean_Fields()
                validDict = djSample.validate()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")
                            
                if prgArgs.upload and validStatus:
                    #_StdProcess.append("Sample")
                    #djCmpd.std_process += ";Sample"
                    djSample.save()
                    outNumbers['Updated Samples'] += 1    
                #------------------------------------------------------------
            

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
