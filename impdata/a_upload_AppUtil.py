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
# from oraCastDB import oraCastDB
# from zUtils import zData

from apputil.models import ApplicationUser, Dictionary

#-----------------------------------------------------------------------------------
def clean_Char(charValue,empty=True):
#-----------------------------------------------------------------------------------
    if charValue:
        return(charValue)
    else:
        if empty:
            return("")
    return(None)

#-----------------------------------------------------------------------------------
def split_List(strList,sep=";"):
#-----------------------------------------------------------------------------------
    if strList:
        retLst = strList.split(sep)
        for i in range(len(retLst)):
            retLst[i] = retLst[i].strip()
    else:
        retLst = None
    return(retLst)
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
def update_AppUser_xls(XlsFile, XlsSheet=0,upload=False):
#-----------------------------------------------------------------------------------
    rmColumns = ['cU','cN','cI','is_staff']
    if os.path.exists(XlsFile):
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}]")
        dfSheet = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        # df -> lstDict and remove null items 
        lstDict = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='rows')]
        for entry in lstDict:

            # remove additional (non-model) columns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]

            # create model entry
            djUser=ApplicationUser(**entry)
            if upload:
                logger.info(f" -> {djUser} ")
                djUser.save()
            else:
                logger.info(f" >r {djUser} ")

#-----------------------------------------------------------------------------------
def update_Dictionary_xls(XlsFile, XlsSheet=0, upload=False, uploaduser=None, lower=False):
#-----------------------------------------------------------------------------------
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

        # df -> lstDict and remove null items 
        lstDict = [{k:v for k,v in m.items() if pd.notnull(v)} for m in dfSheet.to_dict(orient='rows')]

        # check user
        appuser = ApplicationUser.get(uploaduser)

        for entry in lstDict:
            # remove additional (non-model) columns as per rmColumns
            for rmCol in rmColumns:
                if rmCol in entry:
                    del entry[rmCol]

            # check if instance exists in DB
            try:
                djDict = Dictionary.objects.get(dict_value=entry['dict_value'])
            except:
                djDict = Dictionary()

            # set values in instance
            for e in entry:
                setattr(djDict,e,entry[e])
            
            #djDict=Dictionary(**entry)
            if upload:
                logger.info(f" -> {djDict} as {appuser}")
                djDict.save(user=appuser)
            else:
                logger.info(f" >r {djDict} as {appuser}")
