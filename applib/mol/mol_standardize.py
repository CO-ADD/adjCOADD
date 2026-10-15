import os,sys
import json,csv
import logging
import numpy as np
import pandas as pd
import math

from abc import abstractmethod, ABC
from typing import List, Tuple, Dict, Sequence, Optional

from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import rdinchi
from collections import Counter
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import Descriptors

from molvs import Standardizer
#from molvs.metal import MetalDisconnector
#from molvs.normalize import NORMALIZATIONS, MAX_RESTARTS, Normalizer
#from molvs.tautomer import TAUTOMER_TRANSFORMS, TAUTOMER_SCORES, MAX_TAUTOMERS, TautomerCanonicalizer, TautomerEnumerator
#from molvs.charge import ACID_BASE_PAIRS, CHARGE_CORRECTIONS, Reionizer, Uncharger

#from molvs import validate_smiles
#from molvs import standardize_smiles

import logging
logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------------------------
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
AtomType['Metalloids']  = ['B','Si','As','Te','At']
AtomType['Alkali']      = ['Li','Na','K','Rb','Cs','Fr']
AtomType['AlkaliEarth'] = ['Be','Mg','Ca','Sr','Ba','Ra']
AtomType['Halogen']     = ['F', 'Cl','Br','I']
AtomType['Organic']     = ['C','N','O','P','S','Se']
#AtomType['Sulphur']     = ['S']
#AtomType['Selenium']    = ['Se']
#AtomType['MetallAll']   = AtomType['Metall'] + AtomType['MetallTrans'] + AtomType['MetalLanAct']


# # --------------------------------------------------------------------------------------------
# def is_atomtype(at,atype):
#     if atype in AtomType:
#         aSymbol = at.GetSymbol()
#         return (aSymbol in AtomType[atype])
#     return()

# # --------------------------------------------------------------------------------------------
# def list_atomtype_in_mol(mol,atype,unique=True):
#     if atype in AtomType:
#         if mol:
#             alst = []
#             for atom in mol.GetAtoms():
#                 atSym = atom.GetSymbol()
#                 if atSym in AtomType[atype]:
#                     alst.append(atSym)
#             if unique:
#                 alst = list(set(alst))
#             return(alst)
#     return()


# --------------------------------------------------------------------------------------------
class MolStandardizer(ABC):
# --------------------------------------------------------------------------------------------
    """
    Interface for the standardization of molecules, given as SMILES strings.
    """
# --------------------------------------------------------------------------------------------

    def __call__(self, molecules: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
        return self.run(molecules)

    def run(self, molecules: Sequence[str]) -> Sequence[str]:
        """
        Standardizes a sequence of molecules, using run_single
        """
        single_results = [self.run_single(m) for m in molecules]
        return single_results

    @abstractmethod
    def run_single(self, molecule: str) -> str:
        """
        Standardizes a single  molecule.
        """

    def df_apply(self,s,colSMI='SMILES',colSTD='SMI_STD',colVal='VALID'):
        """
        DataFrame.apply
        """
        s[colSTD], s[colVal] = self.run_single(s[colSMI])
        return(s)

# --------------------------------------------------------------------------------------------
class SmiStandardizer_molvs(MolStandardizer):
# --------------------------------------------------------------------------------------------
    """
    Using the MolVS Standardizer

        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = self.disconnect_metals(mol)
        mol = self.normalize(mol)
        mol = self.reionize(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    Args:
        MolStandardizer (_type_): _description_
    """
    def __init__(self,
                 validate: bool = False, uncharge: bool = True, normalize: bool = False, desalt: bool = True,
                 chemdb = None, param_dir: str = None, 
                 ) -> None:
        
        super().__init__()
        self.validate = validate
        self.normalize = normalize
        self.desalt = desalt
        self.uncharge = uncharge
        if param_dir is None:
            param_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'Data')
        if uncharge:
            self.Uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        if normalize:
            self.Normalizer = MolNormalizer(norm_type='chNorm', chemdb=chemdb, param_dir=param_dir)
        if desalt:
            self.Desalter = MolDesalt(chemdb=chemdb, paramdir=param_dir)
            
        self.molvStandardizer = Standardizer()

        
    def run_single_mol(self, smi: Sequence[str]):
        try:
            _mol = Chem.MolFromSmiles(smi)
            _valid = 1
        except:
            _mol = None
            _valid = 0
        if _mol is not None:
            _mol = update_mol_valences(_mol)
            _mol = remove_sgroups_from_mol(_mol)
            _mol = kekulize_mol(_mol)
            _mol = self.molvStandardizer.standardize(_mol)
            _mol = remove_hs_from_mol(_mol)

            if self.desalt:
                _mol= self.Desalter.run_single(_mol)['mol']
            if self.normalize:
                _mol= self.Normalizer.run_single(_mol)
            if self.uncharge:
                _mol = self.Uncharger.uncharge(_mol)
                _mol.UpdatePropertyCache(strict=False)

            if self.validate:
                charge = Chem.GetFormalCharge(_mol)
                if charge > 0:
                    print(f"[Charge]: +{charge} ")
                elif charge < 0:
                    print(f"[Charge]: {charge} ")

            #return(Chem.MolToSmiles(_mol,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False))
            return(_mol,_valid)
        return(None,0)    

    def run_single(self, smi: Sequence[str]) -> str:
        _smi = None
        _mol,_valid = self.run_single_mol(smi)
        if _mol is not None:
            try:
                _smi=Chem.MolToSmiles(_mol,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
            except:
                _smi=Chem.MolToSmiles(_mol,isomericSmiles=True,allHsExplicit=False)
        return(_smi,_valid)

# --------------------------------------------------------------------------------------------
class SmiStandardizer_DB(MolStandardizer):
# --------------------------------------------------------------------------------------------
    """
    Using the MolVS Standardizer

        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
        mol = self.disconnect_metals(mol)
        mol = self.normalize(mol)
        mol = self.reionize(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    Args:
        MolStandardizer (_type_): _description_
    """
    def __init__(self,
                 validate: bool = False, uncharge: bool = True, normalize: bool = False, desalt: bool = True,

                 chemdb = None, param_dir: str = None, 
                 ) -> None:
        
        super().__init__()
        self.validate = validate
        self.normalize = normalize
        self.desalt = desalt
        self.uncharge = uncharge
        if param_dir is None:
            param_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'Data')
        if uncharge:
            self.Uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        if normalize:
            self.Normalizer = MolNormalizer(norm_type='chNorm', chemdb=chemdb, param_dir=param_dir)
        if desalt:
            self.Desalter = MolDesalt(chemdb=chemdb, paramdir=param_dir)
            
        self.molvStandardizer = Standardizer()

    #------------------------------------------------------------------    
    def run_single_mol(self, sMol: Sequence[str]):
    #------------------------------------------------------------------    

        _moldict  = {}
        _saltdict = {}
        _iondict  = {}
        _solvdict = {}

    
        # try:
        #     _mol = Chem.MolFromSmiles(smi)
        # except:
        #     _mol = None
        _mol = sMol
        if _mol is not None:
            _mol = update_mol_valences(_mol)
            _mol = remove_sgroups_from_mol(_mol)
            _mol = kekulize_mol(_mol)

            #_lst_atomtype = list_atomtype_in_mol(_mol)
            try:
                _mol = self.molvStandardizer.standardize(_mol)
                _mol = remove_hs_from_mol(_mol)
                
                if self.desalt:
                    _moldict,_saltdict,_iondict,_solvdict = self.Desalter.run_single(_mol)
                    if _moldict['nfrag'] > 0:
                        _mol = _moldict['mol']
                    else:
                        _mol = None
            except:
                _mol = None

        if _mol is not None:
            if self.normalize:
                _mol= self.Normalizer.run_single(_mol)
            if self.uncharge:
                _mol = self.Uncharger.uncharge(_mol)
                _mol.UpdatePropertyCache(strict=False)

            if self.validate:
                _moldict['charge'] = Chem.GetFormalCharge(_mol)
                # if charge > 0:
                #     print(f"[Charge]: +{charge} ")
                # elif charge < 0:
                #     print(f"[Charge]: {charge} ")

            # Return main Frag Info
            _moldict['mol'] = _mol
            _moldict['valid'] = 1
            _moldict['mw'] = Descriptors.MolWt(_mol)

            # Return Salt/Ion/Solvent - MW
            _moldict['mw_extra'] = 0
            _moldict['smiles_extra'] = ""
            for k,d in _saltdict.items():
                _moldict['mw_extra'] += d['n'] * d['mw']
                if isinstance(d['n'],int):
                    for n in range(d['n']):
                        _moldict['smiles_extra'] += f".{d['smiles']}"
            for k,d in _iondict.items():
                _moldict['mw_extra'] += d['n'] * d['mw']
                if isinstance(d['n'],int):
                    for n in range(d['n']):
                        _moldict['smiles_extra'] += f".{d['smiles']}"
            for k,d in _solvdict.items() :
                _moldict['mw_extra'] += d['n'] * d['mw']
                if isinstance(d['n'],int):
                    for n in range(d['n']):
                        _moldict['smiles_extra'] += f".{d['smiles']}"

            try:
                _moldict['smi']=Chem.MolToSmiles(_moldict['mol'],kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
            except:
                _moldict['smi']=Chem.MolToSmiles(_moldict['mol'],isomericSmiles=True,allHsExplicit=False)

        else:
            _moldict['valid'] = 0
            _moldict['mol'] = None
        return(_moldict,_saltdict,_iondict,_solvdict)


    #------------------------------------------------------------------    
    def run_single(self, smi: Sequence[str]) -> str:
    #------------------------------------------------------------------    

        try:
            _mol = Chem.MolFromSmiles(smi)
        except:
            _mol = None

        _moldict,_saltdict,_iondict,_solvdict = self.run_single_mol(_mol)

        # _smi = None
        # if _moldict['valid']>0:
        #     try:
        #         _smi=Chem.MolToSmiles(_moldict['mol'],kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
        #     except:
        #         _smi=Chem.MolToSmiles(_moldict['mol'],isomericSmiles=True,allHsExplicit=False)
        # else:
        #     _moldict['mol'] = None 

        # if _smi:
        #     _moldict['smi'] = _smi
        return(_moldict,_saltdict,_iondict,_solvdict)
  
# --------------------------------------------------------------------------------------------
class SmiStandardizer(MolStandardizer):
# --------------------------------------------------------------------------------------------
    
    def __init__(self, 
                 validate=False, uncharge=True, normalize=True, desalt=True
                 ) -> None:
        
        super().__init__()
        self.validate = validate
        self.normalize = normalize
        self.desalt = desalt
        self.uncharge = uncharge
        if uncharge:
            self.Uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
        if normalize:
            self.Normalizer = MolNormalizer(normType='chNorm')
        if desalt:
            self.Desalter = MolDesalt()
 
    #-----------------------------------------------------------------------------
    def run_single(self, smi: Sequence[str]) -> str:
        _valid = ""
        
        try:
            _mol = Chem.MolFromSmiles(smi)
            rvalid = 1
        except:
            _mol = None
            rsmi = None
            rvalid = 0
            
        if _mol is not None:
            _mol = update_mol_valences(smol)
            _mol = remove_sgroups_from_mol(smol)
            _mol = kekulize_mol(smol)
            _mol = remove_hs_from_mol(smol)
            if self.desalt:
                molDict= self.Desalter.run_single(_mol)
                smol = molDict['mol']
                #print(molDict)
            if self.normalize:
                smol= self.Normalizer.run_single(_mol)
            if self.uncharge:
                smol = self.Uncharger.uncharge(_mol)

                smol.UpdatePropertyCache(strict=False)
            if self.validate:
                charge = Chem.GetFormalCharge(_mol)
                if charge > 0:
                    print(f"[Charge]: +{charge} ")
                elif charge < 0:
                    print(f"[Charge]: {charge} ")

            rsmi = Chem.MolToSmiles(smol,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
            
        return(rsmi,rvalid)    

#-----------------------------------------------------------------------------
class MolNormalizer(MolStandardizer):
#-----------------------------------------------------------------------------

    def __init__(self,
                 norm_type: str = 'chNorm',
                 chemdb = None, param_dir: str = None, param_file: str = "zChem_Reactions.csv"
                 ) -> None:
        
        super().__init__()    
        self.norm_type = norm_type
        self.alkoxide_pattern = Chem.MolFromSmarts('[Li,Na,K;+0]-[#7,#8;+0]')
        self.normalizer_params = rdMolStandardize.CleanupParameters()
 
        self.param_file = param_file
        self.param_filedir = param_dir
        self.chemdb = chemdb
        
        _normalization_transforms = ""

        if self.chemdb:
            lstReaction = self.chemdb.get_Reaction_List(reaction_type=self.norm_type)
            for n in lstReaction:
                _normalization_transforms = _normalization_transforms + n['smarts'] + '\n'
            logger.info(f"[MolNormalizer] {len(lstReaction)} from ChemDB {self.norm_type}")

        elif self.param_file and self.param_filedir:
            _param_file = os.path.join(self.param_filedir,self.param_file)
            if os.path.isfile(_param_file):
                with open(_param_file) as inFile:
                    csvFile = csv.DictReader(inFile)
                    csvFile.fieldnames = [name.lower() for name in csvFile.fieldnames]
                    for row in csvFile:
                        if row['reaction_type'] == self.norm_type:
                            _normalization_transforms = _normalization_transforms + row['smarts'] + '\n'
                    logger.info(f"[MolNormalizer] {len(row)} from {self.param_file} {self.norm_type}")
        
        self.Normalizer = rdMolStandardize.NormalizerFromData(_normalization_transforms,self.normalizer_params)
        
    def run_single(self,m):
        Chem.FastFindRings(m)
        if m.HasSubstructMatch(self.alkoxide_pattern):
            m = Chem.RWMol(m)
            for match in m.GetSubstructMatches(self.alkoxide_pattern):
                m.RemoveBond(match[0], match[1])
                m.GetAtomWithIdx(match[0]).SetFormalCharge(1)
                m.GetAtomWithIdx(match[1]).SetFormalCharge(-1)

        res = self.Normalizer.normalize(m) 
        return(res)

#-----------------------------------------------------------------------------
class MolDesalt(MolStandardizer):
#-----------------------------------------------------------------------------

    def __init__(self, return_salts = True,
                 chemdb = None, paramdir: str = None, paramfile: str = "zChem_Salt.csv"
                 ) -> None:
        
        super().__init__()    

        self.return_salts = return_salts
        self.param_file = paramfile
        self.param_filedir = paramdir
        self.chemdb = chemdb
        
        self.data = { 'Salt': [], 'Ion': [], 'Solvent': [],}
        
        if self.chemdb:
            for _salt_type in self.data:
                lstSalt = self.chemdb.get_Salt_DictList(_salt_type)
                for row in lstSalt:
                    row['smol'] = Chem.MolFromSmiles(row['smiles'])
                    row['mw'] = Chem.Descriptors.MolWt(row['smol'])
                    self.data[_salt_type].append(row) 
                logger.info(f"[MolDesalt] {_salt_type} {len(self.data[_salt_type])} from ChemDB {chemdb}")

        elif self.param_file and self.param_filedir:
            _param_file = os.path.join(self.param_filedir,self.param_file)
            if os.path.isfile(_param_file):
                with open(_param_file) as inFile:
                    csvFile = csv.DictReader(inFile)
                    csvFile.fieldnames = [name.lower() for name in csvFile.fieldnames]
                    for row in csvFile:
                        if row['salt_type'] in self.data:
                            try:
                                row['smol'] = Chem.MolFromSmiles(row['smiles'])
                                row['mw'] = Chem.Descriptors.MolWt(row['smol'])
                                self.data[row['salt_type']].append(row) 
                            except:
                                logger.warning(f"[MolDesalt] Wrong {row['salt_id']} : {row['smiles']}")
                    logger.info(f"[MolDesalt] {len(row)} from {self.param_file} ")

    def run_single(self,m):
        #desaltDict = get_Fragments_mol(m,self.data,maxDuplicates=1,keepSaltFrag=self.return_salts)
        #print(desaltDict)
        return(get_Fragments_mol(m,self.data,maxDuplicates=1,keepSaltFrag=self.return_salts))
       
#-----------------------------------------------------------------------------
def update_mol_valences(m):
#-----------------------------------------------------------------------------
    m = Chem.Mol(m)
    m.UpdatePropertyCache(strict=False)
    return m

#-----------------------------------------------------------------------------
def remove_sgroups_from_mol(m):
#-----------------------------------------------------------------------------
    # removes all Sgroups
    Chem.ClearMolSubstanceGroups(m)
    return m
#-----------------------------------------------------------------------------
def kekulize_mol(m):
#-----------------------------------------------------------------------------
    Chem.Kekulize(m)
    return m

#-----------------------------------------------------------------------------
def gen_isomers_mol(m, max_isomers=20, embedding=False):
#-----------------------------------------------------------------------------
    
    if m:
        opts = StereoEnumerationOptions(tryEmbedding=embedding,unique=True,maxIsomers=max_isomers,rand=0xf00d)
        nIso = GetStereoisomerCount(m,options=opts)
        #print(f"nIso: {nIso}")
        if nIso <= max_isomers:
            isomers = tuple(EnumerateStereoisomers(m, options=opts))
            return(isomers)
        else:
            return([m])
    return([])

#-----------------------------------------------------------------------------
def gen_isomers_smi(smi, as_smi=True, max_isomers=20, embedding=False):
#-----------------------------------------------------------------------------
    _isosmi = []
    if smi: 
        _mol = Chem.MolFromSmiles(smi)
        _isomers = gen_isomers_mol(_mol, max_isomers=max_isomers, embedding=embedding)
        for _iso in _isomers:
            try:
                _smi=Chem.MolToSmiles(_iso,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
                _isosmi.append(_smi)
            except:
                _smi=Chem.MolToSmiles(_iso,isomericSmiles=True,allHsExplicit=False)
                _isosmi.append(_smi)
    return(_isosmi) 

    #Chem.MolToSmiles(smol,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False)
#-----------------------------------------------------------------------------
def get_tautomers_mol(m, as_smi=True, max_isomers=10):
    te = rdMolStandardize.TautomerCanonicalizer()
    tauts = te.Enumerate(m)
    if as_smi:
        return([Chem.MolToSmiles(m,kekuleSmiles=True,isomericSmiles=True,allHsExplicit=False) for m in tauts])
    return tauts
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def remove_hs_from_mol(m):
#-----------------------------------------------------------------------------
    """ removes most Hs

    Hs that are preserved by the RDKit's Chem.RemoveHs() will not
    be removed.

    Additional exceptions:
    - Hs with a wedged/dashed bond to them
    - Hs bonded to atoms with tetrahedral stereochemistry set
    - Hs bonded to atoms that have three (or more) ring bonds that are not simply protonated
    - Hs bonded to atoms in a non-default valence state that are not simply protonated 


    For the above, the definition of "simply protonated" is an atom with charge = +1 and
    a valence that is one higher than the default.

    """
    # we need ring info, so be sure it's there (this won't do anything if the rings
    # have already been found)
    Chem.FastFindRings(m)
    if m.NeedsUpdatePropertyCache():
        m.UpdatePropertyCache(strict=False)
    SENTINEL = 100
    for atom in m.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree(
        ) == 1 and not atom.GetIsotope():
            nbr = atom.GetNeighbors()[0]
            bnd = atom.GetBonds()[0]
            preserve = False
            if bnd.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH) or \
                    (bnd.HasProp("_MolFileBondStereo") and bnd.GetUnsignedProp("_MolFileBondStereo") in (1, 6)):
                preserve = True
            else:
                is_protonated = nbr.GetFormalCharge() == 1 and \
                    nbr.GetExplicitValence() == \
                    Chem.GetPeriodicTable().GetDefaultValence(nbr.GetAtomicNum())+1
                if nbr.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                                          Chem.ChiralType.CHI_TETRAHEDRAL_CW):
                    preserve = True
                elif not is_protonated:
                    if nbr.GetExplicitValence() > Chem.GetPeriodicTable(
                    ).GetDefaultValence(nbr.GetAtomicNum()):
                        preserve = True
                    else:
                        ringBonds = [
                            b for b in nbr.GetBonds()
                            if m.GetRingInfo().NumBondRings(b.GetIdx())
                        ]
                        if len(ringBonds) >= 3:
                            preserve = True
            if preserve:
                # we're safe picking an arbitrary high value since you can't do this in a mol block:
                atom.SetIsotope(SENTINEL)

    res = Chem.RemoveHs(m, sanitize=False)
    for atom in res.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetIsotope() == SENTINEL:
            atom.SetIsotope(0)
    return res

#-----------------------------------------------------------------------------
def normalize_mol(m,normalizer):
#-----------------------------------------------------------------------------
    _alkoxide_pattern = Chem.MolFromSmarts('[Li,Na,K;+0]-[#7,#8;+0]')
    
    Chem.FastFindRings(m)
    if m.HasSubstructMatch(_alkoxide_pattern):
        m = Chem.RWMol(m)
        for match in m.GetSubstructMatches(_alkoxide_pattern):
            m.RemoveBond(match[0], match[1])
            m.GetAtomWithIdx(match[0]).SetFormalCharge(1)
            m.GetAtomWithIdx(match[1]).SetFormalCharge(-1)
    res = normalizer.normalize(m)
    return res

#-----------------------------------------------------------------------------
def uncharge_mol(m):
#-----------------------------------------------------------------------------
    """
    >>> def uncharge_smiles(smi): return Chem.MolToSmiles(uncharge_mol(Chem.MolFromSmiles(smi)))
    >>> uncharge_smiles('[NH3+]CCC')
    'CCCN'
    >>> uncharge_smiles('[NH3+]CCC[O-]')
    'NCCCO'
    >>> uncharge_smiles('C[N+](C)(C)CCC[O-]')
    'C[N+](C)(C)CCC[O-]'
    >>> uncharge_smiles('CC[NH+](C)C.[Cl-]')
    'CCN(C)C.Cl'
    >>> uncharge_smiles('CC(=O)[O-]')
    'CC(=O)O'
    >>> uncharge_smiles('CC(=O)[O-].[Na+]')
    'CC(=O)[O-].[Na+]'
    >>> uncharge_smiles('[NH3+]CC(=O)[O-].[Na+]')
    'NCC(=O)[O-].[Na+]'
    >>> uncharge_smiles('CC(=O)[O-].C[NH+](C)C')
    'CC(=O)O.CN(C)C'

    Alcohols are protonated before acids:

    >>> uncharge_smiles('[O-]C([N+](C)C)CC(=O)[O-]')
    'C[N+](C)C(O)CC(=O)[O-]'

    And the neutralization is done in a canonical order, so atom ordering of the input
    structure isn't important:

    >>> uncharge_smiles('C[N+](C)(C)CC([O-])CC[O-]')
    'C[N+](C)(C)CC([O-])CCO'
    >>> uncharge_smiles('C[N+](C)(C)CC(CC[O-])[O-]')
    'C[N+](C)(C)CC([O-])CCO'

    """
    uncharger = rdMolStandardize.Uncharger(canonicalOrder=True)
    res = uncharger.uncharge(m)
    res.UpdatePropertyCache(strict=False)
    return res

#-----------------------------------------------------------------------------
def n_Fragments(m):
#-----------------------------------------------------------------------------
    inputFrags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    return(len(inputFrags))

#-----------------------------------------------------------------------------
def find_Fragments(SaltList,sFrag):
#-----------------------------------------------------------------------------
    xDict = {}
    for salt in SaltList:
        for frag in sFrag:
            if 'Salt' not in frag:
                if frag['mol'].HasSubstructMatch(salt['smol']) \
                and salt['smol'].HasSubstructMatch(frag['mol']):
                    if salt['salt_id'] in xDict:
                        #xDict[salt['salt_id']] += 1
                        xDict[salt['salt_id']]['n'] += 1
                    else:
                        xDict[salt['salt_id']] = {'n':1,'smiles':salt['smiles'],'mf':salt['mf'],'mw':salt['mw']}
                    frag['Salt'] = salt['salt_id']
    return(xDict)

def n_Salt_Fragments(sDict):
    _nSalt = 0
    if len(sDict) >0:
        for k,v in sDict.items():
            _nSalt += v['n']
    return(_nSalt)


#-----------------------------------------------------------------------------
def get_Fragments_mol(m,StandardizeList,maxDuplicates=0,keepSaltFrag=False):
#-----------------------------------------------------------------------------

    _subSMatch_param = Chem.SubstructMatchParameters()
    _subSMatch_param.useChirality = True
    _subSMatch_param.useEnhancedStereo = True

    fragList = []
    inputFrags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    for frag in inputFrags:
        fdict = {}
        fdict['mol'] = Chem.RemoveHs(frag, sanitize=False)
        fdict['mol'].UpdatePropertyCache(strict=False)
        Chem.SetAromaticity(fdict['mol'])
        fragList.append(fdict)
    nFrag = len(fragList)

    SaltDict = find_Fragments(StandardizeList['Salt'],fragList)
    IonDict = find_Fragments(StandardizeList['Ion'],fragList)
    SolventDict = find_Fragments(StandardizeList['Solvent'],fragList)
    
    # Check if just Salt/Ion/Solvent
    #print(SaltDict)
    molFrag = nFrag - (n_Salt_Fragments(SaltDict)+n_Salt_Fragments(IonDict)+n_Salt_Fragments(SolventDict))

    MolDict = {}
    if molFrag > 0:
        if nFrag == 1:
            # Single Molecule
            MolDict = fragList[0]
            MolDict['nfrag'] = 1
        else:
            # Multiple Molecules - Find List of non Salt/Ion/Solvent Molecules
            parentList = []
            mwList = []
            smiList = []
            for frag in fragList:
                if 'Salt' not in frag and 'Solvent' not in frag and 'Ion' not in frag:
                    parentList.append(frag['mol'])
            nParentFrag = len(parentList)

            # All Molecules are Salt/Ion/Solvent - Leave largest Salt/Ion/Solvent by nAtoms
            if keepSaltFrag:
                if nParentFrag == 0:
                    #maxMW = 0
                    maxNAtom = 0
                    largFrag = -1
                    nn = -1
                    for frag in fragList:
                        nn += 1
                        xNAtom = frag['mol'].GetNumAtoms()
                        if xNAtom > maxNAtom:
                            maxNAtom = xNAtom
                            largFrag = nn
                        #xMW = Chem.Descriptors.MolWt(frag['mol'])
                        #if xMW > maxMW:
                        #    maxMW = xMW
                        #    largFrag = nn
                    # Reset largest Salt/Ion/Solvent as Parent
                    if largFrag > -1:
                        if 'Salt' in fragList[largFrag]:
                            SaltDict.pop(fragList[largFrag]['Salt'])
                            fragList[largFrag].pop('Salt')
                        if 'Ion' in fragList[largFrag]:
                            IonDict.pop(fragList[largFrag]['Ion'])
                            fragList[largFrag].pop('Ion')
                        if 'Solvent' in fragList[largFrag]:
                            SolventDict.pop(fragList[largFrag]['Solvent'])
                            fragList[largFrag].pop('Solvent')
                    # Reset Parent List
                    parentList = []
                    for frag in fragList:
                        if 'Salt' not in frag and 'Solvent' not in frag and 'Ion' not in frag:
                            parentList.append(frag['mol'])
                    nParentFrag = len(parentList)

            MolDict = {}
            # Single Molecule after Salt/Ion/Solvent
            if nParentFrag == 1:
                MolDict['mol'] = parentList[0]

            # Empty Parent List   
            elif nParentFrag == 0:
                MolDict['mol'] = m
                xFrag = nFrag
                SaltDict = {}
                IonDict = {}
                SolventDict = {}
    
            # Combine multiple Fragment into Parent
            else:
                # Check duplicates
                dupList = []
                for n1 in range(nParentFrag-1):
                    for n2 in range(n1+1,nParentFrag):
                        if parentList[n1].HasSubstructMatch(parentList[n2],_subSMatch_param) \
                        and parentList[n2].HasSubstructMatch(parentList[n1],_subSMatch_param):
                            dupList.append(n2)
                dSalt = (nParentFrag-len(dupList)) / nParentFrag
                MolDict['nDup'] = len(dupList)
                if len(dupList)>0 and len(dupList) <= maxDuplicates:
                    # Remove duplicates
                    #print(dSalt," - ",dupList)
                    dupList.sort(reverse=True)                 
                    for d in dupList:
                        parentList.pop(d)
                    nParentFrag = len(parentList)

                    for s,d in SaltDict.items():
                        SaltDict[s]['n'] = round(d['n']*dSalt,2)
                    for s,d in IonDict.items():
                        IonDict[s]['n'] = round(d['n']*dSalt,2)
                    for s,d in SolventDict.items():
                        SolventDict[s]['n'] = round(d['n']*dSalt,2)

                # Combine multiple molecules     
                if nParentFrag > 1:
                    mix = parentList[0]
                    for sng in parentList[1:]:
                        mix = Chem.CombineMols(mix, sng)
                    MolDict['mol'] = mix
                else:
                    MolDict['mol'] = parentList[0]

            MolDict['nfrag'] = nParentFrag
    else:
        MolDict['nfrag'] = 0

    # MolDict['Salt'] = SaltDict
    # MolDict['Ion'] = IonDict
    # MolDict['Solvent'] = SolventDict

    return(MolDict, SaltDict, IonDict, SolventDict)

