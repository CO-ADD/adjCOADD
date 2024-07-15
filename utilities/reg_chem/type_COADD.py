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



AtomType = {}
AtomType['MetallTrans'] = [
        'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
        'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
        'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn']
AtomType['MetalLanAct'] = [
        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Ac','Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']    
AtomType['Metall']      = ['Al', 'Ga','Ge', 'In','Sn','Sb', 'Tl','Pb','Bi','Po']
AtomType['Alkali']      = ['Li','Na','K','Rb','Cs','Fr']
AtomType['AlkaliEarth'] = ['Be','Mg','Ca','Sr','Ba','Ra']
AtomType['Metalloids']  = ['B','Si','As','Te','At']
AtomType['Halogen']     = ['F', 'Cl','Br','I']
AtomType['Organic']     = ['C','N','O','P','S','Se']

# ==========================================================================
def is_atomtype(at,atype):
# ==========================================================================
    if atype in AtomType:
        aSymbol = at.GetSymbol()
        return (aSymbol in AtomType[atype])
    return()

# ==========================================================================
def list_atomtype_in_mf(mf,atype,unique=True):
# ==========================================================================
    if atype in AtomType:
        if mf:
            alst = []
            for m in AtomType[atype]:
                if len(m)>1:
                    if m in mf:
                        alst.append(m)
            if unique:
                alst = list(set(alst))
            return(alst)
    return()

# ==========================================================================
def list_atomtype_in_mol(mol,atype,unique=True):
# ==========================================================================
    if atype in AtomType:
        if mol:
            alst = []
            for atom in mol.GetAtoms():
                atSym = atom.GetSymbol()
                if atSym in AtomType[atype]:
                    alst.append(atSym)
            if unique:
                alst = list(set(alst))
            return(alst)
    return()

# ==========================================================================
def list_mftype(mf,unique=True):
# ==========================================================================
    mfType = {}
    for atype in AtomType:
        for qatm in AtomType[atype]:
            if qatm in mf:
                mfType[atype] = 1
    return(mfType)

# ==========================================================================
def list_moltype(mol,unique=True):
# ==========================================================================
    molType = {}
    if mol:
        for atom in mol.GetAtoms():
            atSym = atom.GetSymbol()
            for atype in AtomType:
                if atSym in AtomType[atype]:
                    molType[atype] = 1
    return(molType)

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

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "Compound" :


        #MolStd = SmiStandardizer_DB(chemdb=Chem_Salt) 

        print("--> COADD_Compound ---------------------------------------------------------")
        qryCmpd = COADD_Compound.objects.all()
        print(f"[CO-ADD Compound] {len(qryCmpd)}")
        print("-------------------------------------------------------------------------")
        OutFile = f"regChem_COADD_{logTime:%Y%m%d_%H%M%S}.xlsx"

        outNumbers = {'Proc':0,'New Compounds':0,'Updated Compounds':0, 'New Samples': 0, 'Upload Samples': 0}

        for djCmpd in tqdm(qryCmpd):
            outNumbers['Proc'] += 1
            updated_sample = False
            if djCmpd.reg_smiles:
                try:
                    _mol = Chem.MolFromSmiles(djCmpd.reg_smiles)
                    _valid = 1
                except:
                    _mol = None
                    _valid = 0
                if _valid> 0:
                    djCmpd.std_structure_type = list(list_moltype(_mol))
                    updated_sample = True

            elif djCmpd.reg_mf:
                djCmpd.std_structure_type = list(list_mftype(djCmpd.reg_mf))
                updated_sample = True

            validStatus = True
            djCmpd.clean_Fields()
            validDict = djCmpd.validate()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('Warning',k,validDict[k],'-')


            if validStatus and prgArgs.upload and updated_sample:
                outNumbers['Updated Compounds'] += 1
                djCmpd.save()

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
