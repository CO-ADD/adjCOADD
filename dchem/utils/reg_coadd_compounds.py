import os
import pandas as pd
from tqdm import tqdm

from dsample.models import Compound_Batch
from dchem.models import Chem_Salt, Chem_Structure
from applib.mol.smi import get_Structure_Type_Smiles, get_MF_Smiles, SaltDict_to_SaltCode

import logging
logger = logging.getLogger('django')

#------------------------------------------------------------------------------
def standardize_coadd(djCmpd, MolStd=None, outdict={}, upload=False, overwrite=False):
#------------------------------------------------------------------------------
    
    if MolStd is not None:
        
        if not djCmpd.std_status or djCmpd.std_status == 'Invalid' or overwrite:
            _IsMet = 0
            validStatus = False
            updatedStatus = False
            
            if djCmpd.reg_smiles or djCmpd.reg_mf:

                _MolType,_Metal,_IsMet = get_Structure_Type_Smiles(djCmpd.reg_smiles,djCmpd.reg_mf)
                djCmpd.std_structure_type = _MolType
                djCmpd.std_metal = _Metal
                updatedStatus = True
                validStatus = True

            # Excluded from SmiStandardizer - as molvs breaks any metal bonds
            # metal specific Standardizer is required, including OpenSmiles syntax for Metalcomplex
            if _IsMet==1:
                djCmpd.std_status = 'Metal'
                outdict['Metal Compounds'] = outdict.get("Metal Compounds", 0) + 1

                validStatus = True
                updatedStatus = True

            # Non Metal complex structures
            elif djCmpd.reg_smiles:
                outdict['To Standard'] = outdict.get("To Standard", 0) + 1
                _moldict, _saltdict, _iondict, _solvdict = MolStd.run_single(djCmpd.reg_smiles)

                if _moldict['valid'] > 0:
                    djCmpd.std_status = 'Valid'
                    if _moldict['nfrag'] > 1:
                        djCmpd.std_status = 'Mixture'
                        outdict['Mixture'] = outdict.get("Mixture", 0) + 1

                    djCmpd.std_process = "Std"

                    djCmpd.std_smiles = _moldict['smi']
                    djCmpd.std_mw = _moldict['mw']
                    djCmpd.std_nfrag = _moldict['nfrag']
                    djCmpd.std_salt = SaltDict_to_SaltCode(_saltdict)
                    djCmpd.std_ion = SaltDict_to_SaltCode(_iondict)
                    djCmpd.std_solvent = SaltDict_to_SaltCode(_solvdict)
                    djCmpd.std_smiles_extra = _moldict['smiles_extra']
                    djCmpd.std_mw_extra = _moldict['mw_extra']
                    djCmpd.std_mf = get_MF_Smiles(_moldict['smi']+_moldict['smiles_extra'])
                    validStatus = True
                    updatedStatus = True
                else:
                    outdict['Std Failed'] = outdict.get("Std Failed", 0) + 1
                    djCmpd.std_status = 'Invalid'
                    djCmpd.std_process = "Std"
                    validStatus = True
                    updatedStatus = True
            else:
                djCmpd.std_status = 'Empty'
                djCmpd.std_process = "Std"
                outdict['Empty'] = outdict.get("Empty", 0) + 1
                updatedStatus = True
                validStatus = True

            djCmpd.set_defaults_model()
            validDict = djCmpd.validate_fields()

            if validDict:
                validStatus = False
                for k in validDict:
                    logger.warning(f"{k}: {validDict[k]}")

            #print(f" [Cmpd] {djCmpd} {djCmpd.std_status} {validStatus}")

            if validStatus and updatedStatus and upload:
                djCmpd.save()
                outdict['Updated Compounds'] = outdict.get("Updated Compounds", 0) + 1
        else:
            outdict['Already Done'] = outdict.get("Already Done", 0) + 1

#------------------------------------------------------------------------------
def checkissues_coadd(djCmpd, outdict={}, upload=False, overwrite=False):
#------------------------------------------------------------------------------

    if djCmpd.std_status == 'Valid':
        validStatus = True
        updatedStatus = False


        _dmw = djCmpd.reg_mw - djCmpd.std_mw
        _dmf = djCmpd.reg_mf.replace(' ','') != djCmpd.std_mf.replace(' ','')
        
        _has_issue = False
        _issues = []

        if abs(_dmw) > 1:
            _has_issue = True
            if not (djCmpd.std_mw_extra - 1 < _dmw < djCmpd.std_mw_extra + 1):
                outdict['Issues MW'] = outdict.get("Issues MW", 0) + 1
                if djCmpd.reg_mw < 0.5:
                    _issues.append(f"[I] dMW: {_dmw:6.1f} [{djCmpd.std_mw_extra:6.1f}] {djCmpd.reg_mw}")
                if abs(_dmw) >= 2:
                    _issues.append(f"[E] dMW: {_dmw:6.1f} [{djCmpd.std_mw_extra:6.1f}]")
                else:
                    _issues.append(f"[W] dMW: {_dmw:6.1f} [{djCmpd.std_mw_extra:6.1f}]")
            else:
                outdict['Issues Salt'] = outdict.get("Issues Salt", 0) + 1
                _issues.append(f"[w] Salt missing in reg_mw: {_dmw:6.1f}")
        
        if _dmf:
            _has_issue = True
            outdict['Issues MF'] = outdict.get("Issues MF", 0) + 1

            if not djCmpd.reg_mf:
                _issues.append(f"[I] dMF: CxHxNxOx <-> {djCmpd.std_mf}")
            else:
                _issues.append(f"[W] dMF: {djCmpd.reg_mf} <-> {djCmpd.std_mf}")

        if _has_issue:
            djCmpd.std_issues = "; ".join(_issues)
            #logger.warning(f"{djCmpd.compound_id} {_issues}")
            
            djCmpd.set_defaults_model()
            validDict = djCmpd.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    logger.warning(f"{k}: {validDict[k]}")

            if validStatus and upload:
                djCmpd.save()
                outdict['Updated Compounds'] = outdict.get("Updated Compounds", 0) + 1

#------------------------------------------------------------------------------
def updatestd_coadd(djCmpd, outdict={}, upload=False, overwrite=False):
#------------------------------------------------------------------------------

    if djCmpd.std_status == 'Valid':
        validStatus = True
        updatedStatus = False

#------------------------------------------------------------------------------
def regstructure_coadd(djCmpd, outdict={}, upload=False, overwrite=False):
#------------------------------------------------------------------------------

    if djCmpd.std_status == 'Valid':
        validStatus = True
        updatedStatus = False
        new_regchem = True

        if djCmpd.std_nfrag == 1:

            # Check if COAD_Compound has already Structure_ID
            djBatch = Compound_Batch.get(djCmpd.compound_id)
            if djBatch is not None:
                if djBatch.structure_id:
                    new_regchem = False
                    outdict['Has ChemStructures'] = outdict.get("Has ChemStructures", 0) + 1
                    if 'Reg' not in djCmpd.std_process:
                        djCmpd.std_process += '; Reg'
                        updatedStatus = True
                    
            # Process only 'New' Structures 
            if overwrite or new_regchem:
                updated_sample = True
                validStatus = True
                
                #------------------------------------------------------------
                djChem = Chem_Structure.get_bySmiles(djCmpd.std_smiles)
                if djChem is None:
                    djChem = Chem_Structure()
                    djChem.set_molecule(djCmpd.std_smiles)
                    djChem.nfrag = djCmpd.std_nfrag
                    outdict['New ChemStructures'] = outdict.get("New ChemStructures", 0) + 1

                    djChem.set_defaults_model()
                    validDict = djChem.validate_fields()
                    
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            logger.warning(f"{k}: {validDict[k]}")
                            
                    if upload and validStatus:
                        djChem.save()
                        outdict['Updated ChemStructures'] = outdict.get("Updated ChemStructures", 0) + 1
                        if 'Reg' not in djCmpd.std_process:
                            djCmpd.std_process += '; Reg'
                            updatedStatus = True
                else:
                    outdict['Existing ChemStructures'] = outdict.get("Existing ChemStructures", 0) + 1
                    if 'Reg' not in djCmpd.std_process:
                        djCmpd.std_process += '; Reg'
                        updatedStatus = True

                if upload and updatedStatus:
                    djCmpd.save()
                    outdict['Updated Compounds'] = outdict.get("Updated Compounds", 0) + 1

                #------------------------------------------------------------
                djBatch = Compound_Batch.get(djCmpd.compound_id)
                if djBatch is None:
                    djBatch = Compound_Batch()
                    djBatch.compound_id = djCmpd.compound_id
                    djBatch.batch_source = 'COADD'
                    djBatch.batch_code = djCmpd.compound_code

                    outdict['New Batches'] = outdict.get("New Batches", 0) + 1

                djBatch.structure_id = djChem
                djBatch.structure_type = djCmpd.std_structure_type
                _salt_code = []
                if djCmpd.std_salt:
                    _salt_code.append(djCmpd.std_salt)   
                if djCmpd.std_ion:
                    _salt_code.append(djCmpd.std_ion)   
                if djCmpd.std_solvent:
                    _salt_code.append(djCmpd.std_solvent)                            
                djBatch.salt_code = ";".join(_salt_code)
                
                djBatch.smiles_extra = djCmpd.std_smiles_extra
                djBatch.mw_extra = djCmpd.std_mw_extra
                djBatch.full_mw = float(djBatch.mw_extra) + float(djChem.mw)
                djBatch.full_mf = get_MF_Smiles(djCmpd.std_smiles + djBatch.smiles_extra)
                
                djBatch.set_defaults_model()
                validDict = djBatch.validate_fields()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")
                            
                if upload and validStatus:
                    #_StdProcess.append("Sample")
                    #djCmpd.std_process += ";Sample"
                    djBatch.save()
                    outdict['Updated Batches'] = outdict.get("Updated Batches", 0) + 1
                        