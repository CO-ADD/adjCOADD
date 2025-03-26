import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

from dplate.models import MasterPlate, MasterWell, TestPlate, TestWell
#from applib.data.set_fielddata import set_model_arrayfields, set_model_fields, set_model_dicts, set_model_fkeys, set_model_dictarrayfields
from decimal import Decimal


#--------------------------------------------------------------------------------------------------------------
# Add MotherWell to TestWell
#--------------------------------------------------------------------------------------------------------------
def add_mother_to_testwell(djMW,djTW,ClearData=False):
    
    FIELD_LIST = [['cmpbatch_lst','cmpbatch_lst'],
                  ['set_lst','set_lst'],
                  ['test_conc_lst','conc_lst'],
                  ['test_conc_unit_lst','conc_unit_lst'],
                  ['test_conc_type_lst','conc_type_lst'],
                ]
    
    if ClearData:
        djTW.clear_cmpbatch_data()
    
    for f in FIELD_LIST:
        _mw = getattr(djMW, f[0], None) 
        if _mw:
            _tw = getattr(djTW, f[1], None) 
            if _tw:
                _lst = _tw + _mw
            else:
                _lst = _mw
            setattr(djTW,f[1],_lst)
            
    djTW.set_cmpbatch_id()
    
#--------------------------------------------------------------------------------------------------------------
# Add MotherPlate to TestPlate
#--------------------------------------------------------------------------------------------------------------
def add_mother_to_testplate(djMP,djTP,ClearData=False):
    if hasattr(djTP,'wells') and hasattr(djMP,'wells'):
        if djTP.n_wells == djMP.n_wells:
            for w in djMP.wells:
                add_mother_to_testwell(djMP.wells[w],djTP.wells[w],ClearData=ClearData)
        else:
            logger.info(f" [MP->TP] Different n_wells {djMP.n_wells} -> {djTP.n_wells}")     
    else:
        logger.info(f" [MP->TP] Missing well data {djMP.plate_id} -> {djTP.plate_id}")