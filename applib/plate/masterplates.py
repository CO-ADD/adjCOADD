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
from applib.data.set_fielddata import set_arrayFields, set_Fields, set_Dictionaries, set_fkeyFields, set_arrayDictionaries
from decimal import Decimal


# --------------------------------------------------------------------------------
def read_motherplate_prepsheet_xls(xlFile, SheetName='MotherPlates', prefix=None, as_is=False):
# --------------------------------------------------------------------------------
    xlWB = pd.ExcelFile(xlFile)
    xDF = xlWB.parse(SheetName)
    xDF.columns = [c.lower() for c in xDF.columns]

    PropertyList = ['PLATING',
                    'COMPOUND_ID','SET_ID','DILUTION','TEST_CONC','TEST_CONC_UNIT','TEST_SOLVENT_CONC',
                    'COMPOUND2_ID','SET2_ID','DILUTION2','TEST2_CONC','TEST2_CONC_UNIT','TEST2_SOLVENT_CONC']
    lstPl = []
    lstMP = xDF['motherplate_id'].unique()

    grpMP = xDF.groupby(by='motherplate_id')
    for mpid,mpwells in grpMP:

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
            #print(row['motherwell_id'])
            djWell = djMP.get_well(row['motherwell_id'])

            set_Fields(djWell,row,djWell.COPY_FIELDS)
            set_arrayFields(djWell,row,djWell.ARRAY_FIELDS)
            set_Dictionaries(djWell,row,list(djWell.DICTIONARY_FIELDS.keys()))
            djWell.n_cmpbatches = len(djWell.cmpbatch_lst)
            validDict = djWell.check_cmpbatch_id()

            for i in range(len(djWell.test_conc_lst)):
                djWell.test_conc_lst[i] = Decimal(djWell.test_conc_lst[i]).quantize(Decimal("1.0000"))

            if validDict:
                validStatus = False    
                logger.warning(validDict)
                row.update(validDict)
            

        djMP.add_dilutions()
        # for idx,row in mpwells.iterrows():
        #     print
            # validStatus = True
            # if validStatus:

            #     validDict = djWell.check_cmpbatch_id()
            #     if validDict:
            #         validStatus = False    
            #         #print(validDict)
            #         row.update(validDict)

            #     validDict = djWell.check_conc_unit_dictionary()
            #     if validDict:
            #         validStatus = False    
            #         #print(validDict)
            #         row.update(validDict)

            #     djWell.setdefault_fields()
            #     validDict = djWell.validate_fields(exclude=list(arrayFields.keys()))


        #djMP.setdefault_model()
        lstPl.append(djMP)
        logger.info(f"[{djMP.plate_id:25s}] - {djMP.plate_type}  {djMP.n_wells}w  [{_status}]")

    return(lstPl)
