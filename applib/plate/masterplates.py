import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

from dplate.models import MasterPlate, MasterWell
from applib.data.set_fielddata import set_model_arrayfields, set_model_fields, set_model_dicts, set_model_fkeys, set_model_dictarrayfields
from decimal import Decimal


# --------------------------------------------------------------------------------
def read_motherplate_prepsheet_xls(xlFile, SheetName='MotherPlates', prefix=None, as_is=False, **kwargs):
# --------------------------------------------------------------------------------
    xlWB = pd.ExcelFile(xlFile)
    xDF = xlWB.parse(SheetName)
    xDF.columns = [c.lower() for c in xDF.columns]

    lstPl = []
    lstMP = xDF['motherplate_id'].unique()

    grpMP = xDF.groupby(by='motherplate_id')
    for mpid,mpwells in grpMP:
        
        dictPl = {}
        dictPl['plate_id'] = xDF['motherplate_id']
        dictPl['valid_status'] = True
        
        _status = "Exists"
        djMP = MasterPlate.get(mpid, WellData=True, verbose=0)
        if djMP is None:
            #logger.info(f" New Plate {mpid}")
            nWells=384
            djMP = MasterPlate.new(mpid, nWells, PlateType='Mother', WellData=True)
            _status = "New"
        else:
            djMP.load_wells(WellModel=MasterWell)

        djMP.plating = mpwells['plating'].unique()[0]
        djMP.dilution_layout = 'Dilution'

        for idx,row in mpwells.iterrows():
            #print(row['compound_id'])
            if not pd.isna(row['compound_id']) and not pd.isna(row['motherwell_id']):
                # Check if at least compound_id and motherwell_id
                djWell = djMP.get_well(row['motherwell_id'])
                set_model_fields(djWell,row,djWell.COPY_FIELDS)
                set_model_arrayfields(djWell,row,djWell.ARRAY_FIELDS)
                set_model_dicts(djWell,row,list(djWell.DICTIONARY_FIELDS.keys()))
                djWell.n_cmpbatches = len(djWell.cmpbatch_lst)
                validDict = djWell.check_cmpbatch_id()

                if validDict:
                    validStatus = False
                    dictPl['valid_status'] = False    
                    logger.warning(validDict)
                    row.update(validDict)
  
        djMP.add_dilutions()
        djMP.set_defaults_model()

        dictPl['new'] = _status == 'New'
        dictPl['plate'] = djMP
        
        lstPl.append(dictPl)
        logger.info(f"[{djMP.plate_id:25s}] - {djMP.plate_type}  {djMP.n_wells}w  [{_status}]")

    return(lstPl)
