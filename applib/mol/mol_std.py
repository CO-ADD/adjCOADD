import sys, os
import datetime
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
from rdkit import Chem
# from rdkit import Geometry
# from rdkit.Chem import rdinchi
# from rdkit.Chem.MolStandardize import rdMolStandardize
# from rdkit.Chem import rdMolTransforms

# from molvs import Standardizer
# from molvs import validate_smiles
# from molvs import standardize_smiles

# from apputil.models import ApplicationUser, Dictionary, ApplicationLog
# from dchem.models import Chem_Structure
#-----------------------------------------------------------------------------

def Smiles_to_Mol(smi):
    try:
        _mol = Chem.MolFromSmiles(smi)
        _valid = 1
    except:
        _mol = None
        _valid = 0
    return(_mol,_valid)


# ==========================================================================
AtomType = {}
AtomType['MetalTrans'] = [
        'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
        'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
        'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn']
AtomType['MetalLanAct'] = [
        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Ac','Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']    
AtomType['Metal']      = ['Al', 'Ga','Ge', 'In','Sn','Sb', 'Tl','Pb','Bi','Po']
AtomType['Alkali']      = ['Li','Na','K','Rb','Cs','Fr']
AtomType['AlkaliEarth'] = ['Be','Mg','Ca','Sr','Ba','Ra']
AtomType['Halogen']     = ['F', 'Cl','Br','I']
AtomType['Metalloids']  = ['B','Si','As','Te','At']
AtomType['Organic']     = ['C','N','O','P','S','Se']

AtomType_Order = ['MetalTrans','MetalLanAct','Metal','Alkali','AlkaliEarth','Halogen','Metalloids','Organic']
Metal_SMI = ['MetalTrans','MetalLanAct','Metal']

#-----------------------------------------------------------------------------
def list_mftype(mf,unique=True):
#-----------------------------------------------------------------------------
    qrymf = mf
    mfType = {}
    metalLst = {}
    metalSMI = 0
    for atype in AtomType_Order:
        for qatm in AtomType[atype]:
            if qatm in qrymf:
                mfType[atype] = 1
                if 'Metal' in atype:
                    metalLst[qatm] = 1
                qrymf = qrymf.replace(qatm,"")
                if atype in Metal_SMI:
                    metalSMI = 1
    return(mfType,metalLst,metalSMI)

#-----------------------------------------------------------------------------
def list_moltype(mol,unique=True):
#-----------------------------------------------------------------------------
    molType = {}
    metalLst = {}
    metalSMI = 0
    if mol:
        for atom in mol.GetAtoms():
            atSym = atom.GetSymbol()
            for atype in AtomType:
                if atSym in AtomType[atype]:
                    molType[atype] = 1
                    if 'Metal' in atype:
                        metalLst[atSym] = 1
                    if atype in Metal_SMI:
                        metalSMI = 1
    return(molType,metalLst,metalSMI)

#-----------------------------------------------------------------------------
def get_Structure_Type_Smiles(qSmi,qMF,sep=';'):
#-----------------------------------------------------------------------------
    _valid = 0
    retMolType = None
    retMetal = None
    IsMet = 0
    if qSmi:
        try:
            _mol = Chem.MolFromSmiles(qSmi)
            _valid = 1
        except:
            _mol = None
            _valid = 0

    if _valid>0:
        _MolType,_Metal,IsMet = list_moltype(_mol)
        if len(_MolType)>0:
            _molList = list(_MolType)
            _molList.sort()
            retMolType = sep.join(_molList)
        if len(_Metal)>0:
            _metList = list(_Metal)
            _metList.sort()
            retMetal = sep.join(_metList)
    else:
        if qMF:
            _MolType,_Metal,IsMet = list_mftype(qMF)
            if len(_MolType)>0:
                _molList = list(_MolType)
                _molList.sort()
                retMolType = f"{sep.join(_MolType)}; #mf"
            if len(_Metal)>0:
                _metList = list(_Metal)
                _metList.sort()
                retMetal = f"{sep.join(_metList)}; #mf"
    return(retMolType,retMetal,IsMet)

#-----------------------------------------------------------------------------
def get_Structure_Type(qMol,qMF,sep=';'):
#-----------------------------------------------------------------------------
    _valid = 0
    retMolType = None
    retMetal = None
    IsMet = 0
    if qMol:
        _mol = qMol
        _valid = 1
    else:
        _mol = None
        _valid = 0

    if _valid>0:
        _MolType,_Metal,IsMet = list_moltype(_mol)
        if len(_MolType)>0:
            _molList = list(_MolType)
            _molList.sort()
            retMolType = sep.join(_molList)
        if len(_Metal)>0:
            _metList = list(_Metal)
            _metList.sort()
            retMetal = sep.join(_metList)
    else:
        if qMF:
            _MolType,_Metal,IsMet = list_mftype(qMF)
            if len(_MolType)>0:
                _molList = list(_MolType)
                _molList.sort()
                retMolType = f"{sep.join(_MolType)}; #mf"
            if len(_Metal)>0:
                _metList = list(_Metal)
                _metList.sort()
                retMetal = f"{sep.join(_metList)}; #mf"
    return(retMolType,retMetal,IsMet)

# ==========================================================================

#-----------------------------------------------------------------------------
def SaltDict_to_SaltCode(saltDict,sep=';'):
    # sDict - > key (value); key (value)
    if len(saltDict)>0:
        _lst = []
        for k,d in saltDict.items():
            _lst.append(f"{k} ({d['n']})")
        return(sep.join(_lst))
    else:
        return(None)

#-----------------------------------------------------------------------------
def SaltDictList_to_SaltCode(saltDictList,sep=';'):
    # sDict - > key (value); key (value)
    _lst = []
    for saltDict in saltDictList:
        if len(saltDict)>0:
            for k,d in saltDict.items():
                _lst.append(f"{k} ({d['n']})")
    if len(_lst)>0:
        return(sep.join(_lst))
    else:
        return(None)

def get_MF_Smiles(smi):
    _mol,_valid = Smiles_to_Mol(smi)
    mf = ""
    if _mol:
        mf = Chem.rdMolDescriptors.CalcMolFormula(_mol)
    return(mf)


def if_SimpleName(sName):
    numbers = sum(c.isdigit() for c in sName)
    letters = sum(c.isalpha() for c in sName)
    spaces  = sum(c.isspace() for c in sName)
    hyphens = sum(c=='-' for c in sName)
    others  = len(sName) - numbers - letters - spaces - hyphens
    
    if hyphens < 3 or numbers < 3 or others < 3 :
        return(sName)
    else:
        return(None)

# def list_metalatoms(mol,unique=True,):
#     MetalClass   = AtomClass['Metall'] + AtomClass['MetallTrans'] + AtomClass['MetalLanAct']
#     if mol:
#         alst = []
#         for atom in mol.GetAtoms():
#             atSym = atom.GetSymbol()
#             if atSym in MetalClass:
#                 alst.append(atSym)
#         if unique:
#             alst = list(set(alst))
#         return(alst)

# # ==========================================================================
# def get_atomclass_list(mol ,lstClass = AtomClass):
# # ==========================================================================
#     aType = []
#     for l in lstClass:
#         fType=list_atomclass_in_mol(mol,l)
#         if len(fType)>0:
#             aType.append(l)
#     return(aType)
    
