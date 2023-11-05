#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

import django
from django.apps import apps
from django.db.models import ForeignKey, Model

from zUtils import zData

from apputil.models import ApplicationLog
from apputil.utils.data import append_StrList
from dorganism.utils.utils import get_subdir

from apputil.models import ApplicationUser, Dictionary
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock, OrgBatch_Image
from django.utils.text import slugify



#-----------------------------------------------------------------------------------
def rename_OrgName_xls(XlsFile, XlsSheet=0,upload=False,uploaduser=None, lower=False, OutputN=20):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0


    rmColumns = ['chk']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        orgLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(orgLst)
        for org in orgLst:
            if 'strain_code' not in org:
                org['strain_code'] = ''
            rename_OrgName(org['organism_id'],org['organism_name'], org['strain_code'],uploaduser=uploaduser, upload=upload)


#-----------------------------------------------------------------------------------
def rename_OrgName(OrgID, newOrgName, newOrgCode = None, uploaduser=None, upload=False):
#-----------------------------------------------------------------------------------
    djOrg = Organism.get(OrgID)
    djTax = Taxonomy.get(newOrgName,verbose=1)
    if djOrg and djTax:
        oldOrgName = djOrg.organism_name
        djOrg.organism_name = djTax

        djOrg.clean_Fields()
        validDict = djOrg.validate()
        if validDict:
            logger.info(f" XX {djOrg} {validDict} ")
        else:
            # --- Upload ---------------------------------------------------------
            if upload:
                logger.info(f" -> {djOrg} OrgName => {newOrgName} ")
                ApplicationLog.add('Update','rename_OrgName','Info',uploaduser,OrgID,f'Update Organism_Name:: {newOrgName} (old {oldOrgName})','Completed')
                djOrg.save()
            else:
                logger.info(f" >r {djOrg} OrgName => {newOrgName} ")


#-----------------------------------------------------------------------------------
def get_Models_byForeignKey(fkModel):
#-----------------------------------------------------------------------------------
    all_models = apps.get_models()
    
    objModel = None
    if isinstance(fkModel,Model):
        objModel = Model
    if isinstance(fkModel,str):
        for m in all_models:
            if m.__name__ == fkModel:
                objModel = m
    
    # Find models with Author as a foreign key
    models_with_fk = []
    for m in all_models:
        #print(m.__name__)
        for f in m._meta.fields:
            #print(type(f), f.related_model)
            if isinstance(f, ForeignKey) and f.related_model == objModel:
                models_with_fk.append(m)    

    return(models_with_fk)


#-----------------------------------------------------------------------------------
def rename_OrgID_xls(XlsFile, XlsSheet=0, lower=False, OutputN=20,**kwargs):
#-----------------------------------------------------------------------------------

    nProc = {}
    nProc['Saved'] = 0
    nProc['notClass'] = 0
    nProc['notFound'] = 0


    rmColumns = ['chk']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        if lower:
            #change column name to lowercase 
            mvColumns = {}
            for c in dfSheet.columns:
                mvColumns[c] = c.lower()
            logger.info(mvColumns)
            dfSheet = dfSheet.rename(mvColumns,axis='columns') 

        #orgLst = [{k:v for k,v in m.items()} for m in dfSheet.to_dict(orient='records')]
        # df -> lstDict and remove null items 
        orgLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='records')]
        nTotal = len(orgLst)
        

        #rename_OrgID_allBatches('GN_0696','Bacillus velezensis',None,{'Old Strain Code':'Kp BRf 79136','strain_code':'Ba.ve BRf 79136'},**kwargs)
        #rename_OrgID_allBatches('GN_0697','Bacillus safensis',None,{'Old Strain Code':'Kp BRf 79136','strain_code':'Ba.sa BRf 79197'},**kwargs)
        #rename_OrgID_allBatches('GN_1045','Enterococcus faecalis',None,{'Old Strain Code':'Kp BRf 79136','strain_code':'En.fas PK 49'},**kwargs)

        #rename_OrgID_sngBatches('GN_0952_02','Acinetobacter baumannii',None,{},**kwargs)
        #rename_OrgID_sngBatches('GN_0981_01','Serratia marcescens',None,{},**kwargs)
        rename_OrgID_sngBatches('GN_1149_02','Escherichia coli',None,{},**kwargs)

        # for org in orgLst:
        #     if 'strain_code' not in org:
        #         org['strain_code'] = ''
        #     rename_OrgID(org['organism_id'],org['organism_name'], {'strain_code':org['strain_code']}, **kwargs)



#-----------------------------------------------------------------------------------
def rename_OrgID_allBatches(oldOrgID, newOrgName, newOrgID=None, updateDict = None, delete=True, **kwargs):
#-----------------------------------------------------------------------------------

    djOrg = Organism.get(oldOrgID)
    newTax = Taxonomy.get(newOrgName,verbose=1)
    if djOrg and newTax:
        oldOrgName = djOrg.organism_name
        djOrg.organism_name = newTax
        djOrg.organism_id = newOrgID
        djOrg.strain_ids = append_StrList(djOrg.strain_ids,f"[Old ID: {oldOrgID}]")
        if updateDict:
            for k,v in updateDict.items():
                setattr(djOrg,k,v)

        djOrg.save()
        newOrgID = djOrg.organism_id

        # Test
        # newOrgID = 'GP_0317'
        # djOrg = Organism.get(newOrgID)
        # djOrg.strain_ids = append_StrList(djOrg.strain_ids,f"[prev: {newOrgID}]")
        # Test END

        _obj = f'Organism  ({oldOrgID})'
        _desc = f'Update Organism_id: {newOrgID} (old {oldOrgID})'
        ApplicationLog.add('Update','rename_OrgID_allBatches','Info',kwargs['uploaduser'],_obj,_desc,'Completed')

        logger.info(f"[* Rename Organism] {oldOrgID} => {newOrgID} [{djOrg.strain_ids}]")

        # Find all Models with FK to Organism
        org_fk_models = get_Models_byForeignKey('Organism')
        for fk in org_fk_models:
            if fk.__name__ == 'Organism_Batch':
                qryBatchID = Organism_Batch.objects.filter(organism_id=oldOrgID)
                for b in qryBatchID:
                    rename_OrgBatchID(b, None, djOrg, None, **kwargs)
            else:
                qryInst = fk.objects.filter(organism_id=oldOrgID)
                for q in qryInst:
                    q.organism_id = djOrg
                    q.save()
                    _obj = f"{fk.__name__} ({q.pk})"
                    ApplicationLog.add('Update','rename_OrgID_allBatches','Info',kwargs['uploaduser'],_obj,_desc,'Completed')
                    logger.info(f"[ +- Rename {fk.__name__}] {str(q)}")

        # Change old Organism with new ID in strain_ids, and mark as deleted 
        oldOrg = Organism.get(oldOrgID)
        oldOrg.strain_ids = append_StrList(djOrg.strain_ids,f"[New ID: {newOrgID}]")
        oldOrg.save()
        if delete:
            _obj = f'Organism  ({oldOrgID})'
            _desc = f'Marked as Deleted Organism_id: {oldOrgID}'
            ApplicationLog.add('Deleted','rename_OrgID_allBatches','Info',kwargs['uploaduser'],_obj,_desc,'Completed')
            oldOrg.delete()

#-----------------------------------------------------------------------------------
def rename_OrgID_sngBatches(oldOrgBatchID, newOrgName, newOrgID=None, updateDict = None, **kwargs):
#-----------------------------------------------------------------------------------

    
    djBatch = Organism_Batch.get(oldOrgBatchID)
    djOrg = Organism.get(djBatch.organism_id.organism_id)
    oldOrgID = djOrg.organism_id

    newTax = Taxonomy.get(newOrgName,verbose=1)
    if djOrg and newTax and djBatch:
        oldOrgName = djOrg.organism_name
        djOrg.organism_name = newTax
        djOrg.organism_id = newOrgID
        djOrg.strain_ids = append_StrList(djOrg.strain_ids,f"[Split ID: {oldOrgID}]")
        if updateDict:
            for k,v in updateDict.items():
                setattr(djOrg,k,v)

        djOrg.save()
        newOrgID = djOrg.organism_id

        _obj = f'Organism  ({oldOrgID})'
        _desc = f'New Organism_id: {newOrgID} (from {oldOrgID})'
        ApplicationLog.add('Update','rename_OrgID_sngBatches','Info',kwargs['uploaduser'],_obj,_desc,'Completed')

        logger.info(f"[* New Organism] {oldOrgID} => {newOrgID} [{djOrg.strain_ids}]")

        # # Find all Models with FK to Organism
        rename_OrgBatchID(djBatch, None, djOrg, None, delete=False, **kwargs)

        # # Change old Organism with new ID in strain_ids
        oldOrg = Organism.get(oldOrgID)
        oldOrg.strain_ids = append_StrList(djOrg.strain_ids,f"[Split ID: {newOrgID}]")
        oldOrg.save()

#-----------------------------------------------------------------------------------
def rename_OrgBatchID(OrgBatch, newOrgBatchID=None, newOrg=None, updateDict = None, delete=True, **kwargs):
#-----------------------------------------------------------------------------------
    oldOrgBatchID = OrgBatch.orgbatch_id
    oldOrgID = OrgBatch.organism_id.organism_id
    newOrgID = newOrg.organism_id
        
    if not newOrgBatchID:
        if newOrgID:
            newOrgBatchID = oldOrgBatchID.replace(oldOrgID,newOrgID)

    if OrgBatch:
        OrgBatch.organism_id = newOrg
        OrgBatch.batch_notes = append_StrList(OrgBatch.batch_notes,f"[Old ID: {oldOrgBatchID}]")
        OrgBatch.orgbatch_id = newOrgBatchID

        _desc = f'Update OrgBatch_id: {newOrgBatchID} (old {oldOrgBatchID})'
        _obj = f'OrgBatch ({oldOrgBatchID})'
        OrgBatch.save()

        ApplicationLog.add('Update','rename_OrgBatchID','Info',kwargs['uploaduser'],_obj,_desc,'Completed')
        logger.info(f"[ +- Rename Organism_Batch] {oldOrgBatchID} => {newOrgBatchID} ")

        # Find all Models with FK to OrgBatch
        batch_fk_models = get_Models_byForeignKey('Organism_Batch')
        for fk in batch_fk_models:                
            qryInst = fk.objects.filter(orgbatch_id=oldOrgBatchID)
            for q in qryInst:
                logger.info(f"[  +- Change {fk.__name__}] {str(q)} with {OrgBatch}")

                if fk.__name__ == 'OrgBatch_Image':
                    # Rename orgbatch_ID specific values
                    _oldI = q.image_name
                    _old = q.image_file
                    q.image_name = _oldI.replace(oldOrgBatchID,newOrgBatchID)
                    q.image_file = f"images/orgbatch/{get_subdir(q.image_name)}/{q.image_name}"
                    logger.warning(f"[  +- Rename {fk.__name__} --> Rename Image Files {_old} -> {q.image_file}")
                elif fk.__name__ == 'Genome_Sequence':
                    # Rename orgbatch_ID specific values
                    _old = q.seq_name
                    q.seq_name = _old.replace(oldOrgBatchID,newOrgBatchID)
                    q.source_code = q.source_code.replace(oldOrgBatchID,newOrgBatchID)
                    q.source_link = q.source_link.replace(oldOrgBatchID,newOrgBatchID)
                    logger.warning(f"[  +- Rename {fk.__name__} --> Rename Sequence Folder/Files {_old} -> {q.seq_name}")
                q.orgbatch_id = OrgBatch
                q.save()
                _obj = f"{fk.__name__} ({q.pk})"
                ApplicationLog.add('Update','rename_OrgBatchID','Info',kwargs['uploaduser'],_obj,_desc,'Completed')

        # Change old Batch with new ID in batch_notes, and mark as deleted 
        oldBatch = OrgBatch.get(oldOrgBatchID)
        oldBatch.batch_notes = append_StrList(oldBatch.batch_notes,f"[New ID: {newOrgBatchID}]")
        oldBatch.save()
        if delete:
            _obj = f'OrgBatch ({oldOrgBatchID})'
            _desc = f'Marked as Deleted OrgBatch_id: {oldOrgBatchID}'
            ApplicationLog.add('Deleted','rename_OrgBatchID','Info',kwargs['uploaduser'],_obj,_desc,'Completed')
            oldBatch.delete()
 
