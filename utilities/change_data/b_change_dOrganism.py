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
def rename_OrgID(oldOrgID, newOrgName, newOrgID=None, updateDict = None, **kwargs):
#-----------------------------------------------------------------------------------

    djOrg = Organism.get(oldOrgID)
    newTax = Taxonomy.get(newOrgName,verbose=1)
    if djOrg and newTax:
        oldOrgName = djOrg.organism_name
        djOrg.organism_name = newTax
        djOrg.organism_id = newOrgID
        if updateDict:
            for k,v in updateDict.items():
                setattr(djOrg,k,v)

        #djOrg.save()
        #newOrgID = djOrg.organism_id

        #Test
        newOrgID = 'GP_0317'
        djOrg = Organism.get(newOrgID)

        print(f"[* Rename Organism] {oldOrgID} => {newOrgID} ")

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
                    print(f"[ +- Rename {fk.__name__}] {str(q)}")

#-----------------------------------------------------------------------------------
def rename_OrgBatchID(oldBatch, newOrgBatchID=None, newOrg=None, updateDict = None, **kwargs):
#-----------------------------------------------------------------------------------
    oldOrgBatchID = oldBatch.orgbatch_id
    oldOrgID = oldBatch.organism_id.organism_id
    newOrgID = newOrg.organism_id

    if not newOrgBatchID:
        if newOrgID:
            newOrgBatchID = oldOrgBatchID.replace(oldOrgID,newOrgID)

    if oldBatch:

        print(f"[ +- Rename Organism_Batch] {oldOrgBatchID} => {newOrgBatchID} ")

        batch_fk_models = get_Models_byForeignKey('Organism_Batch')
        for fk in batch_fk_models:                
            qryInst = fk.objects.filter(orgbatch_id=oldOrgBatchID)
            for q in qryInst:
                print(f"[  +- Rename {fk.__name__}] {str(q)}")
                if fk.__name__ == 'OrgBatch_Image':
                    print(f"[  +- Rename {fk.__name__}] --> Rename Image Files")
                elif fk.__name__ == 'Genome_Sequence':
                    print(f"[  +- Rename {fk.__name__}] --> Rename Sequence Folder/Files {q.seq_name}")

