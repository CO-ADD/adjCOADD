import sys, os
import datetime
import numpy as np
import pandas as pd

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import rdinchi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdMolTransforms

from molvs import Standardizer
from molvs import validate_smiles
from molvs import standardize_smiles

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from dchem.models import Chem_Structure
#-----------------------------------------------------------------------------

AtomClass = {}
AtomClass['MetallTrans'] = [
        'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
        'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
        'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn']
AtomClass['MetalLanAct'] = [
        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
        'Ac','Th','Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']    
AtomClass['Metall']      = ['Al', 'Ga','Ge', 'In','Sn','Sb', 'Tl','Pb','Bi','Po']
AtomClass['Metalloids']  = ['B','Si','As','Te','At']
AtomClass['Alkali']      = ['Li','Na','K','Rb','Cs','Fr']
AtomClass['AlkaliEarth'] = ['Be','Mg','Ca','Sr','Ba','Ra']
AtomClass['Halogen']     = ['F', 'Cl','Br','I']
AtomClass['Organic']     = ['C','N','O','P','S','Se']
#AtomClass['MetallAll']   = AtomClass['Metall'] + AtomClass['MetallTrans'] + AtomClass['MetalLanAct']

# ==========================================================================
def is_atomclass(atm,aClass, lstClass = AtomClass):
# ==========================================================================
    if aClass in lstClass:
        aSymbol = atm.GetSymbol()
        return (aSymbol in lstClass[aClass])
    return()

# ==========================================================================
def list_atomclass_in_mf(mf, aClass, lstClass=AtomClass,unique=True, ):
# ==========================================================================
    if aClass in lstClass:
        if mf:
            alst = []
            for m in lstClass[aClass]:
                if m in mf:
                    alst.append(m)
            if unique:
                alst = list(set(alst))
            return(alst)
    return()

# ==========================================================================
def list_atomclass_in_mol(mol, aClass, unique=True, lstClass = AtomClass):
# ==========================================================================
    if aClass in lstClass:
        if mol:
            alst = []
            for atom in mol.GetAtoms():
                atSym = atom.GetSymbol()
                if atSym in lstClass[aClass]:
                    alst.append(atSym)
            if unique:
                alst = list(set(alst))
            return(alst)
    return()

def list_metalatoms(mol,unique=True,):
    MetalClass   = AtomClass['Metall'] + AtomClass['MetallTrans'] + AtomClass['MetalLanAct']
    if mol:
        alst = []
        for atom in mol.GetAtoms():
            atSym = atom.GetSymbol()
            if atSym in MetalClass:
                alst.append(atSym)
        if unique:
            alst = list(set(alst))
        return(alst)

# ==========================================================================
def get_atomclass_list(mol ,lstClass = AtomClass):
# ==========================================================================
    aType = []
    for l in lstClass:
        fType=list_atomclass_in_mol(mol,l)
        if len(fType)>0:
            aType.append(l)
    return(aType)
    
