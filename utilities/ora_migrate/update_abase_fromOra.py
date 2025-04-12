#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

import django
#from djCOADD import djOrgDB
from oraCastDB.oraCastDB import openCastDB
from oraABase.oraABase import openABase, get_CompoundBatch


# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ABase"
logDir = "log"
logFileName = os.path.join(logDir,f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

if not os.path.isdir(logDir):
    os.mkdir(logDir)

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    django.setup()

    logging.getLogger().addHandler(logging.FileHandler(logFileName,mode='w'))

    from apputil.models import Dictionary,ApplicationUser
    from applib.data.set_fielddata import set_model_arrayfields, set_model_fields, set_model_dicts, set_model_fkeys, set_model_dictarrayfields
    from dplate.models import Labware, TestPlate, TestWell,MasterPlate, MasterWell
    from dsample.models import Convert_ProjectID, Convert_CompoundID
    from dscreen.models import Screen_Run
    from dorganism.models import Organism_Batch
    from dcell.models import Cell_Batch
    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID
    from update_utils import convert_castdb_compoundid_from_ora
    from applib.mol.mol_std import get_Structure_Type, get_MF_Smiles, SaltDict_to_SaltCode, Smiles_to_Mol, SaltDictList_to_SaltCode
    from dsample.models import ABase_Compound,ABase_Compound_Batch,Compound_Batch
    from dchem.models import Chem_Structure,Chem_Salt

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # For Testing
    # D000092208 (AntiBio_R009) - DR compound_id is null
    # E00095668 (OPXP_R25) - synMIC set_id/conc_type not null
    # HC162-09-21 
    # E00092682 - SC 

    # TestWells -------------------------------------------------------------
   # ABase ChemStructure -------------------------------------------------------------
    # if prgArgs.table == "ABaseCompounds" :
    #     OutName = f"[{prgArgs.table}]"
    #     OutDict = []
    #     OutFile = f"UpdateABaseCmpd_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
    #     OutNumbers = {'Processed':0,
    #                   'New Compound':0, 'Uploaded Compound':0, 
    #                   'New Compound':0, 'Uploaded Compound':0,
    #                   'New Compound':0, 'Uploaded Compound':0,
    #                   'Failed':0}

    #     strSQL = "Select ObjdID, ObjsMolFormula, ObjsMolMassValue, ObjsMolFile  from ChemStruct "
    #     if int(prgArgs.test)>0:
    #         strSQL += f" Fetch First {int(prgArgs.test)} Rows Only "

    #     ABaseDB = openABase()

    #     logger.info(f" {OutName} ... ")
    #     strDF = pd.DataFrame(ABaseDB.get_dict_list(strSQL))
    #     nTotal = len(strDF)
    #     logger.info(f" {OutName} {nTotal} ")
    #     # print("--------------------------------------------------------------------")
    #     # print(f"{OutName} {strDF.columns} ")

    #     for idx,row in tqdm(strDF.iterrows(),  total=len(strDF), desc="ABase Compound"):
    #         NewEntry = False
    #         validStatus = True
    #         OutNumbers['Processed'] += 1
    #         djCmp =  ABase_Compound.get(row['objdid'])
    #         if djCmp is None:
    #             NewEntry = True
    #             djCmp = ABase_Compound()
    #             djCmp.compound_id = row['objdid']
    #             OutNumbers['New Compound'] += 1
            
    #         if row['objsmolfile']:
    #             _molblock = row['objsmolfile'].read()
    #             djCmp.reg_molfile = _molblock
    #             djCmp.reg_mw = row['objsmolmassvalue']
    #             djCmp.reg_mf = row['objsmolformula']
                


    #         djCmp.set_defaults_model()
    #         validDict = djCmp.validate_fields()
    #         if validDict:
    #             validStatus = False
    #             for k in validDict:
    #                 logger.warning('Warning',k,validDict[k],'-')
    #             OutDict.append(row)

    #         if validStatus:
    #             if prgArgs.upload:
    #                 if NewEntry or prgArgs.overwrite:
    #                     OutNumbers['Uploaded'] += 1
    #                     djCmp.save(user=prgArgs.appuser)
    #     ABaseDB.close()
        
    # ABase ABaseBatches -------------------------------------------------------------
    if prgArgs.table == "ABaseBatches" :
        OutDict = []
        OutFile = f"ABaseBatch_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,
                      'New ABaseCompounds':0, 'Uploaded ABaseCompounds':0, 
                      'New CmpBatches':0, 'Uploaded CmpBatches':0,
                      'New ABaseBatches':0, 'Uploaded ABaseBatches':0,
                      'Failed':0}
        
        vSQL   ="""
        SELECT A.OBJDNO OBJDNO
                ,A.OBJDID OBJDID
                ,A.OBJDBATCHREF OBJDBATCHREF
                ,B.OBJDNAME CompoundName
                ,A.STDYID Study_ID
                ,A.RGSTANALYSIS Analysis
                ,A.RGSTCOMPOSITION Compoisition
                ,A.RGSTDRUGNAME DrugName
                ,B.LIBRID Library_ID
                ,A.DICTMOLID_SALT Salt_ID
                ,A.RGSTSALTEQUIVS Salt_Equiv
                ,A.DICTMOLID_SOLVATE Solvate_ID
                ,A.RGSTSOLVATEEQUIVS Solvate_Equiv
                ,A.RGSTFULLMOLFORMULA Full_MF
                ,A.RGSTFULLMOLMASSVALUE Full_MW
                ,C.OBJSMOLFORMULA MF
                ,C.OBJSMOLMASSVALUE MW
                ,C.OBJSMOLFILE MOLFILE
                ,B.OBJDORIGINATOR Originator
                ,B.OBJDLABNOTEBOOKNO Lab_Notebook_Number
                ,B.DICTID_OBJDSUPPLIER Supplier
                ,A.RGSTSUPPLIEROBJID Supplier
                ,A.RGSTSUPPLIEROBJID Supplier_CatNo
                ,A.RGSTSUPPLIERBATCHREF Supplier_Batch
                ,A.RGSTDATERECEIVED Date_Received
                ,B.DICTID_OBJDTYPE ObjdType
                ,B.OBJDDATEREGISTERED Registered_Date
                ,B.OBJDQTYINITVALUE Init_Value
                ,B.DICTID_QTYUNIT Init_Value_Unit
                -- ,A.RGSTSTATUS RGSTSTATUS ,A.RGSTOBJSOURCETYPE RGSTOBJSOURCETYPE ,A.DICTTABLECODE_RGSTPURPOSE DICTTABLECODE_RGSTPURPOSE
                -- ,A.DICTID_RGSTPURPOSE DICTID_RGSTPURPOSE ,A.DICTMOLIND_SALT DICTMOLIND_SALT,A.DICTMOLIND_SOLVATE DICTMOLIND_SOLVATE
                -- ,A.RGSTFULLMOLMASSDECPLCSENTRD RGSTFULLMOLMASSDECPLCSENTRD,A.DICTTABLECODE_RGSTCONTAINER DICTTABLECODE_RGSTCONTAINER
                --,A.DICTID_RGSTCONTAINER DICTID_RGSTCONTAINER,A.RGSTCONTAINERID RGSTCONTAINERID,A.RGSTWELLPOS RGSTWELLPOS
                --,A.RGSTPURCHASEORDERNO RGSTPURCHASEORDERNO,A.RGSTSOLVENT RGSTSOLVENT,A.RGSTSUPPLIERBATCHREF RGSTSUPPLIERBATCHREF
                --,A.RGSTDATERECEIVED RGSTDATERECEIVED,A.RGSTREGISTRARUSERID RGSTREGISTRARUSERID,A.DELETEDSTATUS DELETEDSTATUS
                --,A.AUDITCREATEDUSERID AUDITCREATEDUSERID,A.AUDITCREATEDTIME AUDITCREATEDTIME,A.AUDITMODIFIEDUSERID AUDITMODIFIEDUSERID
                --,A.AUDITMODIFIEDTIME AUDITMODIFIEDTIME,B.USGPID USGPID
                --,C.OBJSMOLMASSDECPLCSENTRD PARENT_MASS_DEC_PLACES
        FROM OBJDRGST A
            Left Join OBJDTAIL B on A.OBJDNO = B.OBJDNO
            Left Join CHEMSTRUCT C on  A.OBJDID = C.OBJDID
        """

        ABaseDB = openABase()
        OutName = f"[{prgArgs.table}]"
        cmpData = ABaseDB.get_dict_list(vSQL,upCase=True)
        logger.info(f" {OutName} ABase Batches: {len(cmpData)}") 
        for row in tqdm(cmpData, desc='ABase'):
            OutNumbers['Processed'] += 1
            cbid = f"{row['OBJDID']}_{row['OBJDBATCHREF']}"
            
            # ABase Compound -------------------------------------
            djAbaseCmp = ABase_Compound.get(row['OBJDID'])
            NewCompound = False
            if djAbaseCmp is None:
                NewCompound = True
                djAbaseCmp = ABase_Compound()
                djAbaseCmp.compound_id = row['OBJDID']
                OutNumbers['New ABaseCompounds'] += 1
            
            if row['MOLFILE']:
                _molblock = row['MOLFILE'].read()
                djAbaseCmp.reg_molfile = _molblock
                djAbaseCmp.reg_mw = row['MW']
                djAbaseCmp.reg_mf = row['MF']

            # try:
            #     _aMol = Chem.MolFromMolBlock(_molblock)
            #     Chem.Kekulize(_aMol)
            # except:
            #     _aMol = None

            # if _aMol is not None:
            #     djCmp.reg_smiles = Chem.MolToSmiles(_aMol)
            #     _mf_rdkit = rdMolDescriptors.CalcMolFormula(_aMol)
            #     _mass_rdkit = Descriptors.ExactMolWt(_aMol)
            #     _mw_rdkit = Descriptors.MolWt(_aMol)

            djAbaseCmp.set_defaults_model()
            validDict = djAbaseCmp.validate_model(verbose=0)
                
            if prgArgs.upload:
                if NewCompound or prgArgs.overwrite:
                    OutNumbers['Uploaded ABaseCompounds'] += 1
                    djAbaseCmp.save()


            # CmpBatch -------------------------------------
            djCmpBatch = Compound_Batch.get(cbid)
            NewCmpBatch = False
            if djCmpBatch is None:
                NewCmpBatch = True
                OutNumbers['New CmpBatches'] += 1
                djCmpBatch = Compound_Batch()
                djCmpBatch.cmpbatch_id = cbid
                            
            djCmpBatch.batch_notes = row['DRUGNAME']
            djCmpBatch.batch_id = row['OBJDBATCHREF']
            # djCmpBatch.salt_code
            djCmpBatch.full_mw = row['FULL_MW']
            djCmpBatch.full_mf = row['FULL_MF']
            djCmpBatch.batch_code = f"{row['OBJDID']}:{row['OBJDBATCHREF']}"
            djCmpBatch.batch_source = 'ABASE'

            djCmpBatch.set_defaults_model()
            validDict = djCmpBatch.validate_model(verbose=0)
            
            if prgArgs.upload:
                if NewCmpBatch or prgArgs.overwrite:
                    OutNumbers['Uploaded CmpBatches'] += 1
                    djCmpBatch.save()

            # ABase Compound Batch -------------------------------------
            djAbaseBatch = ABase_Compound_Batch.get(cbid)
            NewABase = False
            if djAbaseBatch is None:
                NewABase = True
                OutNumbers['New ABaseBatches'] += 1
                djAbaseBatch = ABase_Compound_Batch()
                djAbaseBatch.cmpbatch_id = djCmpBatch
                djAbaseBatch.compound_id = djAbaseCmp
        
            #djPrj = Project.get()
            #djAbaseBatch.project_id = djPrj
            djAbaseBatch.library_id = row['LIBRARY_ID']
            
            djAbaseBatch.full_mw = row['FULL_MW']
            djAbaseBatch.full_mf = row['FULL_MF']   
            djAbaseBatch.salt_code = row['SALT_ID']   
            djAbaseBatch.salt_equivalents  = row['SALT_EQUIV']     
            djAbaseBatch.solvate_code = row['SOLVATE_ID']      
            djAbaseBatch.solvate_equivalents = row['SOLVATE_EQUIV']      

            djAbaseBatch.supplier = row['SUPPLIER']        
            djAbaseBatch.supplier_code  = row['SUPPLIER_CATNO']       
            djAbaseBatch.supplier_batch = row['SUPPLIER_BATCH']        
            djAbaseBatch.date_recieved   = row['DATE_RECEIVED']
            
            djAbaseBatch.init_amount = row['INIT_VALUE']   
            djAbaseBatch.init_amount_unit = Dictionary.get(djAbaseBatch.DICTIONARY_FIELDS['init_amount_unit'],row['INIT_VALUE_UNIT'])
            
            
            # Chemist ---------------------------
            USER_RENAME = {
               'X.Chemist': 'orgdb',
               'A.BadilloVega': 'A.Kavanagh',
               'Ciara.Davis':'C.Davis' 
            }
            djUser = ApplicationUser.get(row['ORIGINATOR'])
            if djUser is None:
                for k in USER_RENAME:
                    if row['ORIGINATOR'] == k:
                        djUser = ApplicationUser.get(USER_RENAME[k])           
            if djUser is None:
                logger.error(f" [Chemist] {row['ORIGINATOR']} not found")
            else:
                djAbaseBatch.chemist = djUser

            if row['LAB_NOTEBOOK_NUMBER'] is not None:
                _lab = row['LAB_NOTEBOOK_NUMBER'].split(chr(160))    
                djAbaseBatch.labbook_no = _lab[0]  
                djAbaseBatch.labbook_page = _lab[1]   
                djAbaseBatch.labbook_page_line = _lab[2]  
                 
            djAbaseBatch.set_defaults_model()
            validDict = djAbaseBatch.validate_model(verbose=1)

            if prgArgs.upload:
                if NewABase or prgArgs.overwrite:
                    OutNumbers['Uploaded ABaseBatches'] += 1
                    djAbaseBatch.save()
                
                
        ABaseDB.close()
        logger.info(f"[{prgArgs.table}] {OutNumbers}")
        
   # ABase Structure -------------------------------------------------------------
    if prgArgs.table == "ABaseStructure" :
        OutName = f"[{prgArgs.table}]"

        OutDict = []
        OutFile = f"ABaseStr_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New':0, 'Uploaded':0, 'Failed':0,'Empty':0, 'not STD':0}


        print(f" {OutName} ---------------------------------------------------------")
        logger.info(f" {OutName} ... ")
        if int(prgArgs.test)>0:
            qryStr = ABase_Compound.objects.all()[:int(prgArgs.test)]
            nEntry = int(prgArgs.test)
        else:
            qryStr = ABase_Compound.objects.all()
            nEntry = qryStr.count()
        logger.info(f" {OutName}: {nEntry} ")

        for djMCC in tqdm(qryStr, total=nEntry, desc='[CmpBatcheLsts]'):
            OutNumbers['Processed'] += 1

            if djMCC.reg_molfile:
                _mol = Chem.MolFromMolBlock(djMCC.reg_molfile)
                if _mol:
                    validStatus = True
                    Chem.Kekulize(_mol)
                    _MolType,_Metal,_IsMet = get_Structure_Type(_mol,None)
                    djMCC.structure_type = _MolType
                    djMCC.structure_metal = _Metal

                    _smi = Chem.MolToSmiles(_mol, canonical=True, isomericSmiles=True, kekuleSmiles=True)

                    djChem = Chem_Structure.get_exact(_mol)
                    if djChem is None:
                        djChem = Chem_Structure()
                        djChem.smol = _mol
                        OutNumbers['New'] += 1

                    djChem.set_defaults_model()
                    validDict = djChem.validate_fields()
                    
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            logger.warning(f"{k}: {validDict[k]}")
                            
                    if prgArgs.upload and validStatus:
                        djChem.save()
                        _csid = djChem.structure_id
                        OutNumbers['Uploaded'] += 1 
                        # Reload to get MW
                        djChem.get(_csid)
                        djMCC.structure_id = djChem

                else:
                    _mol = Chem.MolFromMolBlock(djMCC.reg_molfile, sanitize=False)
                    if _mol:
                        OutNumbers['not STD'] += 1
                        _MolType,_Metal,_IsMet = get_Structure_Type(_mol,None)
                        djMCC.structure_type = _MolType
                        djMCC.structure_metal = _Metal
                    else:    
                        OutNumbers['Failed'] += 1

                validStatus = True
                validDict = djMCC.validate_fields()
                djMCC.set_defaults_model()
                validDict = djMCC.validate_fields()
                
                if validDict:
                    validStatus = False
                    for k in validDict:
                        logger.warning(f"{k}: {validDict[k]}")

                if prgArgs.upload and validStatus:
                    djMCC.save()

            else:
                OutNumbers['Empty'] += 1
        logger.info(f"[{prgArgs.table}] {OutNumbers}")
        


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [ABaseCompounds/ABaseBatches/ABaseStructure]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    from zDjango.djUtils import init_django_dir

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)

#==============================================================================
