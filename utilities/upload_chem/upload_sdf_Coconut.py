#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse
from rdkit import Chem 
from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    # Logger ----------------------------------------------------------------
    import logging
    logTime= datetime.datetime.now()
    logName = "Upload_SDF_Coconut"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    #    handlers=[logging.StreamHandler()],
        level=logLevel)
    #-----------------------------------------------------------------------------

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from dchem.models import Chem_Structure, Chem_Salt
    from dchem.utils.mol_std import if_SimpleName
    from dsample.models import Library, Library_Compound

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    LibraryID = 'COCONUT'

    # Table -------------------------------------------------------------
    if prgArgs.table == "Library":
        
        appuser = ApplicationUser.get(prgArgs.appuser)
        djLib = Library.get(LibraryID)
        if djLib:
            # with Chem.ForwardSDMolSupplier(prgArgs.file) as sdSupl:
            #     lines = 0
            #     for mol in sdSupl:
            #         lines += 1

            outNumbers = {'Proc':0,'New Compounds':0,'Upload Compounds':0, 'New Samples': 0, 'Upload Samples': 0, 'Failed': 0}
            outDict = []    
            with Chem.ForwardSDMolSupplier(prgArgs.file) as sdSupl:
                #reader = csv.DictReader(infile)
#                for mol in tqdm(sdSupl,total=lines, desc="Reading Library"):
                for mol in tqdm(sdSupl, desc="Reading SDFile"):
                    outNumbers['Proc'] += 1
                    if mol is not None:
                        new_compound = False
                        row = mol.GetPropsAsDict()
                        compound_code = row['coconut_id']

                        djCmpd = Library_Compound.get(None,compound_code,LibraryID)
                        if not djCmpd:
                            djCmpd = Library_Compound()
                            djCmpd.compound_code = compound_code
                            djCmpd.library_id = djLib
                            djCmpd.compound_name = if_SimpleName(row['name'],250)
                            djCmpd.compound_desc = f"Found in: {row['found_in_databases']}"
                            djCmpd.reg_smiles = row['SMILES']
                            new_compound = True

                        #set_dictFields(djCmpd,row,['compound_name','reg_smiles'])
                        
                        validStatus = True

                        djCmpd.clean_Fields()
                        validDict = djCmpd.validate()
                        if validDict:
                            validStatus = False
                            for k in validDict:
                                logger.warning(f"{k} : {validDict[k]}")
                            outDict.append(row)

                        if validStatus:
                            if prgArgs.upload:
                                if new_compound or prgArgs.overwrite:
                                    outNumbers['Upload Compounds'] += 1
                                    djCmpd.save()
                    else:
                        outNumbers['Failed'] += 1
            logger.info(f"[LibCompounds] :{outNumbers}")

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