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

TypeConf = {
            'None':['Unknown','Synthetic'],
            'Antibiotic':['Small molecule','Synthetic'],
            'Extract':['Extract','Natural product'],
            'Extracts':['Extract','Natural product'],
            'Metal complex of Small Synthetic molecule':['Transmetal complex','Synthetic'],
            'Mixture; Peptide':['Peptide','Synthetic'],
            'Mixture; Peptide; Small Molecule':['Peptide','Synthetic'],
            'Mixture; Small Molecule':['Small molecule','Synthetic'],
            'Mixture; Small Molecule; Peptide':['Peptide','Synthetic'],
            'Nanoparticle':['Nanoparticle','Synthetic'],
            'Natural Product':['Small molecule','Natural product'],
            'Natural Product; Derivative':['Small molecule','Semisynthetic'],
            'Natural Product; Extract':['Extract','Natural product'],
            'Natural Product; Peptide':['Peptide','Natural product'],
            'Organometallic Molecule':['Organometallic molecule','Synthetic'],
            'Peptide':['Peptide','Synthetic'],
            'Peptide; Cyclic':['Peptide cyclic','Synthetic'],
            'Peptide; Cyclic; Amphiphile':['Peptide cyclic','Synthetic'],
            'Peptidemimetic':['Peptidemimetic','Synthetic'],
            'Polymer':['Polymer','Synthetic'],
            'Reference product1':['Small molecule','Synthetic'],
            'Reference product2':['Small molecule','Synthetic'],
            'Reference product3':['Small molecule','Synthetic'],
            'Reference product5':['Small molecule','Synthetic'],
            'Reference product7':['Small molecule','Synthetic'],
            'Semisynthetic':['Small molecule','Semisynthetic'],
            'Small Molecule':['Small molecule','Synthetic'],
            'Small  Molecule':['Small molecule','Synthetic'],
            'Small Molecule; Polar':['Small molecule','Synthetic'],
            'small synthetic molecule':['Small molecule','Synthetic'],
            'Small synthetic molecule':['Small molecule','Synthetic'],
            'Small Synthetic molecule':['Small molecule','Synthetic'],
            'Small Synthetic Molecule':['Small molecule','Synthetic'],
            'Solvent':['Solvent','Synthetic'],
            'Synthetic':['Small molecule','Synthetic'],
            'Synthetic compound':['Small molecule','Synthetic'],
            'Synthetic Product':['Small molecule','Synthetic'],
            'Synthetic Reference product4':['Small molecule','Synthetic'],
            'Transmetal Complex':['Transmetal complex','Synthetic'],
            'Transmetal Complex; Peptide; Cyclic':['Transmetal complex','Synthetic'],
            'Transmetal Complex; Polymer':['Transmetal complex','Synthetic'],
            'Transmetal Salt':['Salt only','Synthetic'],
        }



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
        qryCmpd = COADD_Compound.objects.all()            
        nCmpd = qryCmpd.count()    
        logger.info(f"[CO-ADD Compound] {nCmpd} [All]")
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
            
            ctype = 'None'
            if djCmpd.compound_type:
                _ctype = djCmpd.compound_type

            if djCmpd.ora_compound_type in TypeConf:
                setattr(djCmpd,'compound_type',Dictionary.get(djCmpd.Choice_Dictionary['compound_type'],TypeConf[djCmpd.ora_compound_type][0]))
                setattr(djCmpd,'compound_source',Dictionary.get(djCmpd.Choice_Dictionary['compound_source'],TypeConf[djCmpd.ora_compound_type][1]))

                djCmpd.clean_Fields()
                validDict = djCmpd.validate()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")
                            
                if prgArgs.upload and validStatus:
                    #_StdProcess.append("Sample")
                    #djCmpd.std_process += ";Sample"
                    djCmpd.save()
                    outNumbers['Updated Compounds'] += 1    
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
