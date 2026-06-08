
import sys, os
import datetime
import numpy as np
import pandas as pd
from decimal import Decimal, getcontext
from tqdm import tqdm

import logging

from apputil.models import AuditModel
from dchem.models import Chem_Structure
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
#-----------------------------------------------------------------------------

from apputil.models import Dictionary, ApplicationUser
from applib.external.sql_oracle import Oracle
from dsample.models import Compound_Batch, ABase_Compound_Batch, ABase_Compound, Project


#-----------------------------------------------------------------------------
def openABase():
    db = Oracle()
    db.open("ABase","ABASE","imb-coadd-db.imb.uq.edu.au","1521","coadb")
    return(db)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def get_ABase_Batches(test=0):
#-----------------------------------------------------------------------------
    
    ABASE_SQL   ="""
    SELECT A.OBJDNO OBJDNO
            ,A.OBJDID OBJDID 
            ,A.OBJDBATCHREF 
            ,B.OBJDNAME
            ,A.STDYID Study_ID
            ,A.RGSTANALYSIS Analysis
            ,A.RGSTCOMPOSITION Compoisition
            ,A.RGSTDRUGNAME Drug_Name
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
 
    if test>0:
        ABASE_SQL += f" Fetch First {test} Rows Only "
       
    ABaseDB = openABase()
    print(f"[ABase Batches] ... ")
    _DF = pd.DataFrame(ABaseDB.get_dict_list(ABASE_SQL))
    nTotal = len(_DF)
    print(f"[ABase Batches] {nTotal} ")
    ABaseDB.close()
    return(_DF)

#-----------------------------------------------------------------------------
def get_ABase_RegView(test=0):
#-----------------------------------------------------------------------------

    _SQL = "Select * from ObjdRgst_View "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        _SQL += f" Fetch First {test} Rows Only "

    ABaseDB = openABase()
    print(f"[ABaseRegView] ... ")
    _DF = pd.DataFrame(ABaseDB.get_dict_list(_SQL))
    nTotal = len(_DF)
    print(f"[ABaseRegView] {nTotal} ")
    ABaseDB.close()

    return(_DF)

#-----------------------------------------------------------------------------
def get_ABaseChem_Structure(CompoundID=None, RDKit=False):
#-----------------------------------------------------------------------------
    strSQL = "Select ObjdID, ObjsMolFormula, ObjsMolMassValue, ObjsMolFile  from ChemStruct "
    
    ABaseDB = openABase()
    
    if CompoundID is not None:
        strSQL += f" Where ObjdID = '{CompoundID}' "
          
    _structList = ABaseDB.get_dict_list(strSQL)
    
    if _structList and len(_structList)>0:
        for i in range(len(_structList)):
            if pd.isna(_structList[i]['objsmolfile']):
                _structList[i]['molfile'] = ""
            else:
                _molfile = _structList[i]['objsmolfile'].read()
                _structList[i]['molfile'] = _molfile
            
                if RDKit:
                    try:
                        aMol = Chem.MolFromMolBlock(_molfile)
                        Chem.Kekulize(aMol)
                    except:
                        aMol = None
                    if aMol is not None:
                        #print(ObjdID," - ",aMol.GetNumAtoms())
                        _structList[i]['rdkit_smiles'] = Chem.MolToSmiles(aMol)
                        _structList[i]['rdlit_mf'] = rdMolDescriptors.CalcMolFormula(aMol)
                        _structList[i]['rdkit_mass'] = Descriptors.ExactMolWt(aMol)
                        _structList[i]['rdkit_mw'] = Descriptors.MolWt(aMol)
                

    ABaseDB.close()
    return(_structList)

#-----------------------------------------------------------------------------
def upload_ABase_Batches(regDF,upload=False,overwrite=False):
#-----------------------------------------------------------------------------

    OutNumbers = {'Processed':0,'New CmpBatch':0,'New ABase':0,'New ABase Batch':0,'Uploaded Entries':0}
    
    print(f"[ABaseBatchDict] {regDF.columns.tolist()} ")
    
    # Set global precision to 6 significant digits
    #getcontext().prec = 6
    
    for idx,row in tqdm(regDF.iterrows(), total=regDF.shape[0], desc='Abase Batches Upload'):
    #for idx,row in regDF.iterrows():
        #print(row)
        OutNumbers['Processed'] += 1
        NewEntry = False
        validStatus = True
        oraBatch_id = f"{row['objdid']}:{row['objdbatchref']}"
        
        # ------------------------------------------------------
        djCompound_id = f"{row['objdid']}"
        djBatch_id = f"{row['objdid']}_{row['objdbatchref']}"


        # ABase Compound ----------------------------------------------------------------
        NewEntry = False
        djABaseCmp = ABase_Compound.get(djCompound_id)
        if djABaseCmp  is None:
            NewEntry = True
            OutNumbers['New ABase'] += 1
            djABaseCmp = ABase_Compound()
            djABaseCmp.compound_id = djCompound_id
            
        _structList = get_ABaseChem_Structure(CompoundID=djCompound_id)
        if _structList and len(_structList)>0: 
            _struct = _structList[0]
            djABaseCmp.reg_mf = _struct['objsmolformula']
            
            #print(f" * [{_struct['objsmolmassvalue']}] [{type(_struct['objsmolmassvalue'])}]")
            # print(_struct['objsmolmassvalue'])
            # print(type(_struct['objsmolmassvalue']))
            if pd.isna(_struct['objsmolmassvalue']):
                #print(_struct['objsmolmassvalue'])
                djABaseCmp.reg_mw = 0
            else:
                #print(type(_struct['objsmolmassvalue']))
                #djABaseCmp.reg_mw = round(Decimal(str(_struct['objsmolmassvalue'])),3)
                djABaseCmp.reg_mw = Decimal(str(_struct['objsmolmassvalue']))
                
            #print(f" * [{_struct['objsmolmassvalue']}] [{djABaseCmp.reg_mw}]")
            djABaseCmp.reg_molfile = _struct['molfile']            
            
        # - Save -------------
        djABaseCmp.set_defaults_model()
        validDict = djABaseCmp.validate_fields()
        if validDict:
            validStatus = False
            for k in validDict:
                logger.warning('Warning',k,validDict[k],'-')
        if validStatus:
            if upload:
                if NewEntry or overwrite:
                    djABaseCmp.save()


        # Cmpound Batch ----------------------------------------------------------------
        NewEntry = False
        djCmpBatch = Compound_Batch.get(djBatch_id)
        if djCmpBatch  is None:
            NewEntry = True
            djCmpBatch = Compound_Batch()
            djCmpBatch.cmpbatch_id = djBatch_id
            djCmpBatch.batch_id = row['objdbatchref']
            OutNumbers['New CmpBatch'] += 1

        djCmpBatch.full_mf = row['full_mf']
        djCmpBatch.full_mw = round(Decimal(row['full_mw']),3)
        djCmpBatch.batch_source = 'ABASE'
        djCmpBatch.batch_code = f"{row['objdid']}:{row['objdbatchref']}"
        if 'rgstdrugname' in row:
            djCmpBatch.batch_notes = row['drug_name']

        # - Save -------------
        djCmpBatch.set_defaults_model()
        validDict = djCmpBatch.validate_fields()
        if validDict:
            validStatus = False
            for k in validDict:
                logger.warning('Warning',k,validDict[k],'-')        
        if validStatus:
            if upload:
                if NewEntry or overwrite:
                    OutNumbers['Upload Entries'] += 1
                    djCmpBatch.save()

       # ABase Compound ----------------------------------------------------------------
        NewEntry = False
        djABaseCmpBatch = ABase_Compound_Batch.get(djBatch_id)
        if djABaseCmpBatch  is None:
            NewEntry = True
            OutNumbers['New ABase Batch'] += 1
            djABaseCmpBatch = ABase_Compound_Batch()
            djABaseCmpBatch.cmpbatch_id = djCmpBatch

        # Study ID ---------------------------
        STUDYID_RENAME = {
            'G01_Antibact': '026_Antibiotic',
        }
        djPrj = Project.objects.filter(abase_study_id = row['study_id']).first()
        if djPrj is None:
            for k in STUDYID_RENAME:
                if row['study_id'] == k:
                    djPrj = Project.objects.filter(abase_study_id = STUDYID_RENAME[k]).first()    
        if djPrj:
            djABaseCmpBatch.project_id = djPrj
        else:
            logger.error(f" [Project] {row['study_id']} not found ")
            
        djABaseCmpBatch.library_id = row['library_id']
        djABaseCmpBatch.full_mw = round(Decimal(row['full_mw']),3)
        djABaseCmpBatch.full_mf = row['full_mf']
           
        djABaseCmpBatch.salt_code = row['salt_id']
        if  pd.isna(row['salt_equiv']) :  
            djABaseCmpBatch.salt_equivalents  = 0
        else:
            djABaseCmpBatch.salt_equivalents  = Decimal(row['salt_equiv'])

        djABaseCmpBatch.solvate_code = row['solvate_id']      
        if  pd.isna(row['solvate_equiv']) :  
            djABaseCmpBatch.solvate_equivalents  = 0
        else:
            djABaseCmpBatch.solvate_equivalents  = Decimal(row['solvate_equiv'])

        if djABaseCmpBatch.full_mw > 0 and djABaseCmp.reg_mw > 0:
            djABaseCmp.conv_factor = round(djABaseCmpBatch.full_mw / djABaseCmp.reg_mw,4)
        else:
            djABaseCmp.conv_factor = 0
        #print(f" [{djABaseCmpBatch.full_mw}] [{djABaseCmp.reg_mw}] -> [{djABaseCmp.conv_factor}]")

        djABaseCmpBatch.supplier = row['supplier']        
        djABaseCmpBatch.supplier_code  = row['supplier_catno']       
        djABaseCmpBatch.supplier_batch = row['supplier_batch']        
        djABaseCmpBatch.date_recieved   = row['date_received']
        
        #print(f" {row['init_value']} {row['init_value_unit']} ")
        djABaseCmpBatch.init_amount = row['init_value']
        djUnit = Dictionary.get(djABaseCmpBatch.DICTIONARY_FIELDS['init_amount_unit'],row['init_value_unit'])
        if djUnit:   
            djABaseCmpBatch.init_amount_unit = djUnit
        else:
            djABaseCmpBatch.init_amount_unit = None
            #logger.error(f" [Unit] {row['init_value_unit']} not found ")

        if row['lab_notebook_number'] is not None:
            _lab = str(row['lab_notebook_number']).split(chr(160))
            djABaseCmpBatch.labbook_no = _lab[0]
            if len(_lab) > 1:  
                djABaseCmpBatch.labbook_page = _lab[1]
                if len(_lab) > 2:   
                    djABaseCmpBatch.labbook_page_line = _lab[2]
            
        # Chemist ---------------------------
        USER_RENAME_CHANGE = {
            'X.Chemist': 'orgdb',
            'A.BadilloVega': 'A.Kavanagh',
            'Ciara.Davis':'C.Davis' 
        }
        djUser = ApplicationUser.get(row['originator'])
        if djUser is None:
            for k in USER_RENAME_CHANGE:
                if row['originator'] == k:
                    djUser = ApplicationUser.get(USER_RENAME_CHANGE[k])           
        if djUser is None:
            logger.error(f" [Chemist] {row['originator']} not found")
        else:
            djABaseCmpBatch.chemist = djUser            
        
            
        # - Save -------------
        djABaseCmpBatch.set_defaults_model()
        # validDict = djABaseCmpBatch.validate_fields()
        # if validDict:
        #     validStatus = False
        #     for k in validDict:
        #         logger.warning('Warning',k,validDict[k],'-')
        # if validStatus:
        if upload:
            if NewEntry or overwrite:
                djABaseCmpBatch.save()


    print(f"[ABaseBatchDict] {OutNumbers}")
    
#-----------------------------------------------------------------------------
def get_ABase_Tests(db):
#-----------------------------------------------------------------------------
    selSQL = """
        SELECT * 
        FROM ABASE.Test
        """

#-----------------------------------------------------------------------------
def get_ABase_TestOccation(db,TOC):
#-----------------------------------------------------------------------------
    selSQL = """
        SELECT 
            A.TTORDirNo,
            A.TTORObjSeqNo,
            A.TTORSeqNo,
            A.TTORDirNo_Parent,
            A.TTORObjSeqNo_Parent,
            A.TTORSeqNo_Parent,
            A.StdyId,
            A.ProtId,
            A.ProtVersionNo,
            A.TOccId,
            A.TRTsId,
            A.DictId_Rslt,
            A.DictId_RsltUnit,
            A.TTORStatus,
            A.TTORRsltAlpha,
            B.DictId_Cond,
            B.CondQualifierAlpha,
            B.CondQualifierValue,
            B.DictNumDecPlaces,
            B.DictId_CondUnit,
            A.ObjDId,
            A.ObjDBatchRef,
            D.RgstFullMolMassValue,
            A.OLPtId,
            A.OLPtPlateInd,
            A.TTORWellReference,
            A.CHILD_DictId_Rslt,
            A.CHILD_DictId_RsltUnit,
            A.CHILD_TTORStatus,
            A.CHILD_TTORRsltAlpha,
            C.DictId_Cond,
            C.CondQualifierAlpha,
            C.CondQualifierValue,
            C.DictNumDecPlaces,
            C.DictId_CondUnit,
            A.CHILD_TTORSeqNo,
            A.CHILD_PrPmValSeqNo_Rslt
        FROM 
            ABASE.TOTSODRS_CHILD_VIEW@AbCooper A,
            ABASE.RSLTCONDVAL_VIEW@AbCooper B,
            ABASE.RSLTCONDVAL_VIEW@AbCooper C,
            ABASE.OBJDRGST@AbCooper D
        WHERE
            (A.TOccId=?) AND
            (A.CondGroupNo_Rslt=B.CondGroupNo_Rslt(+) AND A.PrPmNo_Rslt=B.PrPmNo_Rslt(+)) AND
            (A.CHILD_CondGroupNo_Rslt=C.CondGroupNo_Rslt(+) AND A.CHILD_PrPmNo_Rslt=C.PrPmNo_Rslt(+)) AND
            (A.ObjDId = D.ObjDId(+) AND A.ObjDBatchRef = D.ObjDBatchRef(+))
        """
    selSQL += f" AND (A.TOccId={TOC})"
    # unqTestID := TTORDIRNO . "-" . TTOROBJSEQNO;
    # unqCOND := DICTID_COND . "=" . CONDQUALIFIERALPHA . '=' . DICTID_CONDUNIT;
    # 
    # unqResult[#i] := DICTID_RSLT[#i] . "=" . TTORRSLTALPHA[#i] . "=" . DICTID_RSLTUNIT[#i];
    # unqResult[#nVal+#i] := CHILD_DICTID_RSLT[#i] . "=" . CHILD_TTORRSLTALPHA[#i] . "=" . CHILD_DICTID_RSLTUNIT[#i];
    # remove('TTORRSLTALPHA');
    # remove('DICTID_RSLT');
    # remove('DICTID_RSLTUNIT');
    # remove('CHILD_TTORRSLTALPHA');
    # remove('CHILD_DICTID_RSLT');
    # remove('CHILD_DICTID_RSLTUNIT');


#-----------------------------------------------------------------------------
def get_ABase_TestRequest(db,TRS):
#-----------------------------------------------------------------------------
    selSQL = """
        SELECT 
            A.TTORDirNo,
            A.TTORObjSeqNo,
            A.TTORSeqNo,
            A.TTORDirNo_Parent,
            A.TTORObjSeqNo_Parent,
            A.TTORSeqNo_Parent,
            A.StdyId,
            A.ProtId,
            A.ProtVersionNo,
            A.TOccId,
            A.TRTsId,
            A.DictId_Rslt,
            A.DictId_RsltUnit,
            A.TTORStatus,
            A.TTORRsltAlpha,
            B.DictId_Cond,
            B.CondQualifierAlpha,
            B.CondQualifierValue,
            B.DictNumDecPlaces,
            B.DictId_CondUnit,
            A.ObjDId,
            A.ObjDBatchRef,
            D.RgstFullMolMassValue,
            A.OLPtId,
            A.OLPtPlateInd,
            A.TTORWellReference,
            A.CHILD_DictId_Rslt,
            A.CHILD_DictId_RsltUnit,
            A.CHILD_TTORStatus,
            A.CHILD_TTORRsltAlpha,
            C.DictId_Cond,
            C.CondQualifierAlpha,
            C.CondQualifierValue,
            C.DictNumDecPlaces,
            C.DictId_CondUnit,
            A.CHILD_TTORSeqNo,
            A.CHILD_PrPmValSeqNo_Rslt
        FROM 
            ABASE.TOTSODRS_CHILD_VIEW@AbCooper A,
            ABASE.RSLTCONDVAL_VIEW@AbCooper B,
            ABASE.RSLTCONDVAL_VIEW@AbCooper C,
            ABASE.OBJDRGST@AbCooper D
        WHERE
            (A.CondGroupNo_Rslt=B.CondGroupNo_Rslt(+) AND A.PrPmNo_Rslt=B.PrPmNo_Rslt(+)) AND
            (A.CHILD_CondGroupNo_Rslt=C.CondGroupNo_Rslt(+) AND A.CHILD_PrPmNo_Rslt=C.PrPmNo_Rslt(+)) AND
            (A.ObjDId = D.ObjDId(+) AND A.ObjDBatchRef = D.ObjDBatchRef(+))
        """
    selSQL += f" AND (A.TRTsId={TRS})"
    
# #-----------------------------------------------------------------------------
# def upload_ABase_RegView(regDF,upload=False,overwrite=False):
# #-----------------------------------------------------------------------------

#     OutNumbers = {'Processed':0,'New CmpBatch':0,'New ABase':0,'New ABase Batch':0,'Uploaded Entries':0}
    
#     print(f"[ABaseRegDict] {regDF.columns.tolist()} ")
        
#     for idx,row in tqdm(regDF.iterrows(), total=regDF.shape[0], desc='AbaseRegView Upload'):
#     #for idx,row in regDF.iterrows():
#         #print(row)
#         OutNumbers['Processed'] += 1
#         NewEntry = False
#         validStatus = True
#         oraBatch_id = f"{row['objdid']}:{row['objdbatchref']}"
        
#         # ------------------------------------------------------
#         djCompound_id = f"{row['objdid']}"
#         djBatch_id = f"{row['objdid']}_{row['objdbatchref']}"


#         # ABase Compound ----------------------------------------------------------------
#         djABaseCmp = ABase_Compound.get(djCompound_id)
#         if djABaseCmp  is None:
#             NewEntry = True
#             djABaseCmp = ABase_Compound()
#             djABaseCmp.compound_id = djCompound_id
#             _structList = get_ABaseChem_Structure(CompoundID=djCompound_id)
#             if _structList and len(_structList)>0: 
#                 _struct = _structList[0]
#                 djABaseCmp.reg_mf = _struct['objsmolformula']
#                 djABaseCmp.reg_mw = _struct['objsmolmassvalue']
#                 djABaseCmp.reg_molfile = _struct['molfile']            
#             OutNumbers['New ABase'] += 1
        

#         # Cmpound Batch ----------------------------------------------------------------
#         djCmpBatch = Compound_Batch.get(djBatch_id)
#         if djCmpBatch  is None:
#             NewEntry = True
#             djCmpBatch = Compound_Batch()
#             djCmpBatch.cmpbatch_id = djBatch_id
#             djCmpBatch.batch_id = row['objdbatchref']
#             OutNumbers['New CmpBatch'] += 1

#         djCmpBatch.full_mf = row['rgstfullmolformula']
#         djCmpBatch.full_mw = row['rgstfullmolmassvalue']
#         djCmpBatch.batch_source = 'ABASE'
#         djCmpBatch.batch_code = f"{row['objdid']}:{row['objdbatchref']}"
#         if 'rgstdrugname' in row:
#             djCmpBatch.batch_notes = row['rgstdrugname']

#         djCmpBatch.set_defaults_model()
#         validDict = djCmpBatch.validate_fields()
#         if validDict:
#             validStatus = False
#             for k in validDict:
#                 logger.warning('Warning',k,validDict[k],'-')
#             #OutDict.append(row)
#         #print(f" {validStatus} {prgArgs.upload}")
        
#         if validStatus:
#             if upload:
#                 if NewEntry or overwrite:
#                     OutNumbers['Upload Entries'] += 1
#                     djCmpBatch.save()


#        # ABase Compound ----------------------------------------------------------------
#         djABaseCmpBatch = ABase_Compound_Batch.get(djBatch_id)
#         if djABaseCmpBatch  is None:
#             NewEntry = True
#             djABaseCmpBatch = ABase_Compound_Batch()
#             djABaseCmpBatch.cmpbatch_id = djCmpBatch
#             # djABaseCmp.compound_id = djABaseCmp

#             djABaseCmpBatch.library_id = row['library_id']
#             djABaseCmpBatch.project_id = row['study_id']
            
#             djABaseCmpBatch.full_mw = row['rgstfullmolmassvalue']
#             djABaseCmpBatch.full_mf = row['rgstfullmolformula']   
#             djABaseCmpBatch.salt_code = row['dictmolid_salt']   
#             djABaseCmpBatch.salt_equivalents  = row['rgstsaltequivs']     
#             djABaseCmpBatch.solvate_code = row['dictmolid_solvate']      
#             djABaseCmpBatch.solvate_equivalents = row['rgstsolvateequivs']      

#             # djABaseCmp.conv_factor = row['objdbatchref']

#             djABaseCmpBatch.supplier = row['supplier']        
#             djABaseCmpBatch.supplier_code  = row['rgstsupplierobjid']       
#             djABaseCmpBatch.supplier_batch = row['rgstsupplierbatchref']        
#             djABaseCmpBatch.date_recieved   = row['rgstdatereceived']
            
#             djABaseCmpBatch.init_amount = row['objdqtyinitvalue']
#             djUnit = Dictionary.get(djABaseCmpBatch.DICTIONARY_FIELDS['init_amount_unit'],row['dictid_qty_unit'])
#             if djUnit:   
#                 djABaseCmpBatch.init_amount_unit = djUnit

#             if row['objlabnotebookno'] is not None:
#                 _lab = row['objlabnotebookno'].split(chr(160))    
#                 djABaseCmpBatch.labbook_no = _lab[0]  
#                 djABaseCmpBatch.labbook_page = _lab[1]   
#                 djABaseCmpBatch.labbook_page_line = _lab[2]  


#             # Chemist ---------------------------
#             USER_RENAME_CHANGE = {
#                'X.Chemist': 'orgdb',
#                'A.BadilloVega': 'A.Kavanagh',
#                'Ciara.Davis':'C.Davis' 
#             }
#             djUser = ApplicationUser.get(row['originator'])
#             if djUser is None:
#                 for k in USER_RENAME_CHANGE:
#                     if row['originator'] == k:
#                         djUser = ApplicationUser.get(USER_RENAME_CHANGE[k])           
#             if djUser is None:
#                 logger.error(f" [Chemist] {row['originator']} not found")
#             else:
#                 djABaseCmpBatch.chemist = djUser            
            
#             OutNumbers['New ABase Batch'] += 1

#     print(f"[ABaseRegDict] {OutNumbers}")


# class Compound_Batch(AuditModel):
#     """
#     List of Compound Batches 
#     """
# #-------------------------------------------------------------------------------------------------
#     DICTIONARY_FIELDS = {
#         'batch_type':'CmpBatch_Type',
#     }

#     ID_SEQUENCE = 'CmpBatch'
#     ID_PREFIX = 'CB'
#     ID_PAD = 9
    
#     cmpbatch_id = models.CharField(max_length=15, primary_key=True, verbose_name = "CmpBatch ID")

#     batch_id  = models.CharField(default= '00',max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
#     batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

#     batch_source = models.CharField(max_length=25, choices=CMPBATCH_SOURCES, blank=False, verbose_name = "Batch Source")
#     batch_code = models.CharField(max_length=150, blank=True, verbose_name = "Batch Code")
#     previous_ids = models.CharField(max_length=100, blank=True, verbose_name = "Previous IDs")
#     structure_type = models.CharField(max_length=400, blank=True, verbose_name = "Type")
#     structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
#         db_column="structure_id", related_name="%(class)s_structure_id")
#     salt_code = models.CharField(max_length=120, blank=True, verbose_name = "Salts")
#     smiles_extra = models.CharField(max_length=256, blank=True, verbose_name = "Smiles Salt")
#     mw_extra = models.DecimalField(default=0, max_digits=12, decimal_places=3, verbose_name = "MW Salt")
    
#     full_mw = models.FloatField(default=0, blank=True, verbose_name = "Full MW")
#     full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")

# class ABase_Compound(AuditModel):
#     """
#     List of Abase Compounds as per Registration
#     """
# #-------------------------------------------------------------------------------------------------
#     DICTIONARY_FIELDS = {
#     }

#     ID_SEQUENCE = 'ABase_Compound'
#     ID_PREFIX = 'MCC'
#     ID_PAD = 6

#     compound_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Compound ID")
#     compound_code = models.CharField(max_length=50, blank=True, verbose_name = "Code")
#     compound_name = models.CharField(max_length=250, blank=True, verbose_name = "Name")
#     compound_desc = models.CharField(max_length=250, blank=True, verbose_name = "Comment")

#     reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")
#     reg_molfile = models.TextField(max_length=15, blank=True, verbose_name = "Reg Molfile")
#     reg_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Reg MW")
#     reg_mf = models.CharField(max_length=100, blank=True, verbose_name = "Reg MF")
    
#     structure_type = models.CharField(max_length=400, blank=True, verbose_name = "Type")
#     structure_metal = models.CharField(max_length=100, blank=True, verbose_name = "Std Metal")
#     structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
#         db_column="structure_id", related_name="%(class)s_structure_id")


# class ABase_Compound_Batch(AuditModel):

#     DICTIONARY_FIELDS = {
#         'init_amount_unit':'Unit_Amount',
#     }

#    cmpbatch_id = models.OneToOneField(Compound_Batch, primary_key=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
#                                         db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")
#     compound_id = models.ForeignKey(ABase_Compound, null=True, blank=True, verbose_name = "Compound ID", on_delete=models.DO_NOTHING,
#         db_column="compound_id", related_name="%(class)s_compound_id")

#     library_id = models.CharField(max_length=20, blank=True, verbose_name = "Library ID")
#     project_id = models.ForeignKey(Project, null=True, blank=True, verbose_name = "Project ID", on_delete=models.DO_NOTHING,
#             db_column="project_id", related_name="%(class)s_project_id")

#     full_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Full MW")
#     full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")
#     salt_code = models.CharField(max_length=50, blank=True, verbose_name = "Salt Code")
#     salt_equivalents = models.DecimalField(max_digits=7, decimal_places=2, default=0, verbose_name = "Salt Eq")
#     solvate_code = models.CharField(max_length=50, blank=True, verbose_name = "Solvate Code")
#     solvate_equivalents = models.DecimalField(max_digits=7, decimal_places=2, default=0, verbose_name = "Solvate Eq")
#     conv_factor = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Conv Factor")
#     supplier = models.CharField(max_length=50, blank=True, verbose_name = "Supplier")
#     supplier_code = models.CharField(max_length=50, blank=True, verbose_name = "Supplier Code")
#     supplier_batch = models.CharField(max_length=50, blank=True, verbose_name = "Supplier Batch")
#     # supplier_po
#     date_recieved = models.DateField(null=True, blank=True, verbose_name = "Received")
#     init_amount = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Init Amount")
#     init_amount_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Init Amount Unit", on_delete=models.DO_NOTHING,
#         db_column="init_amount_unit", related_name="%(class)s_init_amount_unit")

#     chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
#         db_column="chemist", related_name="%(class)s_chemist")
#     labbook_no = models.CharField(max_length=10, blank=True, verbose_name = "LabBook")
#     labbook_page = models.CharField(max_length=10, blank=True, verbose_name = "LabBook Page")
#     labbook_page_line = models.CharField(max_length=10, blank=True, verbose_name = "LabBook Page Line")
