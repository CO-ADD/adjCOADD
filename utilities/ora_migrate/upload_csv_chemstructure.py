#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from tqdm import tqdm
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
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "D:/Code/zdjCode/adjCOADD"
    # uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    # if prgArgs.database == 'Work':
    #     djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # elif prgArgs.database == 'WorkLinux':
    #     djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"

    # xlFiles = {
    #     'Application': "ApplicationData_v05.xlsx",
    #     'Drug': "DrugData_v04.xlsx",
    #     'MIC': "LMIC_Data_v06.xlsx",
    #     'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    # }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser
    from dchem.models import Chem_Structure
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "UploadOrgDB"
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

   # Table -------------------------------------------------------------
    ExcelFile = "D:/Upload/CastDB/pgData/Chem_Structure.csv"

    choiceTables = ['ChemStructure',
                    ]
    if prgArgs.table in choiceTables:

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

        if prgArgs.table == 'ChemStructure' and prgArgs.file:

            appuser = ApplicationUser.get(prgArgs.appuser)

            with open(prgArgs.file, mode='r') as csv_file:
                lines = len(csv_file.readlines())
            
            outDict = []    
            with open(prgArgs.file, mode='r') as infile:
                reader = csv.DictReader(infile)
                for row in tqdm(reader,total=lines, desc="Reading ChemStructures"):
                    _nfrag = int(row['nfrag'])
                    if  _nfrag< 2:    
                        if Chem_Structure.exists_bySmiles(row['smiles']):
                            o = Chem_Structure.get_bySmiles(row['smiles'])
                            if o.structure_id != row['structure_id']:
                                row['Issue'] = f"Exists"
                                row['IssueValue'] = f"{o.structure_id}; {o.get_smiles()}"
                                outDict.append(row)
                                print(f" Exists: {row['structure_id']} {row['smiles']} -> {o.structure_id} {o.get_smiles()}")
                        else:
                            c = Chem_Structure()
                            c.set_molecule(row['smiles'])
                            c.nfrag = _nfrag
                            c.structure_id = row['structure_id']
                            c.atom_classes = get_atomclass_list(c.smol)
                            lstMetall=list_metalatoms(c.smol)
                            if len(lstMetall)>0:
                                row['Issue'] = f"Metal"
                                row['IssueValue'] = f"{lstMetall}"
                                outDict.append(row)
                                print(f" Metal: {row['structure_id']} {row['smiles']} -> Metal: {lstMetall}") 
                            else:    
                                if prgArgs.upload:
                                    c.save(user=appuser)

                    else:
                        row['Issue'] = f"nFrag"
                        row['IssueValue'] = f"{_nfrag}"
                        outDict.append(row)
                        print(f" nFrag: {row['structure_id']} {row['smiles']} -> nFrag: {_nfrag}") 
 
            outDF = pd.DataFrame(outDict)
            outDF.to_excel('UploadChemStructure_Issues.xlsx')
    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
