import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

from dplate.models import MasterPlate, MasterWell, Labware
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


# --------------------------------------------------------------------------------
def read_barcode_csv(csvFile, add_current_location=True, **kwargs):
# --------------------------------------------------------------------------------

    def apply_current_location(s):
        if s['BARCODE'] == 'NO READ':
            s['CURRENT_PLATE_ID'] = '-'
            s['CURRENT_WELL_ID'] = '-'
            s['ACTION'] = 'EMPTY'
        else:
            _tube = MasterWell.get(None,None,s['BARCODE'])
            if _tube:
                s['CURRENT_PLATE_ID'] = str(_tube.plate_id)
                s['CURRENT_WELL_ID'] = str(_tube.well_id)

                if ((s['CURRENT_PLATE_ID'] == s['PLATE_ID']) and (s['CURRENT_WELL_ID'] == s['WELL_ID'])):
                    s['ACTION'] = 'SAME'
                # elif (s['CURRENT_PLATE_ID'] == s['PLATE_ID']):
                #     s['ACTION'] = 'SWAP'
                else:
                    s['ACTION'] = 'MOVE'
            else:
                s['CURRENT_PLATE_ID'] = '-'
                s['CURRENT_WELL_ID'] = '-'
                s['ACTION'] = 'NEW'
        return(s)

    dfBarcode = pd.read_csv(csvFile)
    dfBarcode.columns = ['PLATE_ID', 'WELL_ID', 'BARCODE']
    if add_current_location:
        dfBarcode = dfBarcode.apply(apply_current_location,axis=1)
    
    return(dfBarcode)             


# --------------------------------------------------------------------------------
def get_dummy_rack(RackID='MP_DUMMY', PlateSize=96, LabewareID='FLUIDX_10mL',upload=False):
# --------------------------------------------------------------------------------
    RACK_TYPE = 'Storage'
    _plate_id = f"{RackID}_{PlateSize}"
    if MasterPlate.exists(_plate_id):
        return(MasterPlate.get(_plate_id))
    else:
        _MP = MasterPlate().new(_plate_id,PlateSize,RACK_TYPE, WellData=False)
        _MP.labware_id = Labware.get(LabewareID)
        if upload:
            _MP.set_defaults_model()
            _MP.save()
            _MP.remove()
        return(_MP)

# --------------------------------------------------------------------------------
def update_barcode_location(csvFile, LabewareID = 'FLUIDX_10mL', upload=False, verbose=0):
# --------------------------------------------------------------------------------
    RACK_SIZE = 96
    RACK_TYPE = 'Storage'

    dfBC = read_barcode_csv(csvFile)

    DummyRack = get_dummy_rack(upload=upload)
    DummyRackID = str(DummyRack.plate_id)

    # Get all Target Rack ----------------------------------------
    TargetRackID = dfBC['PLATE_ID'].unique()[0]
    if MasterPlate.exists(TargetRackID):
        Racks = {TargetRackID:MasterPlate.get(TargetRackID)}
    else:
        Racks = {TargetRackID:MasterPlate().new(TargetRackID,RACK_SIZE,RACK_TYPE, WellData=True)}
        Racks[TargetRackID].labware_id = Labware.get(LabewareID)
        if upload:
            Racks[TargetRackID].save()

    # Get all Source Racks ----------------------------------------
    RackIDs = dfBC['CURRENT_PLATE_ID'].unique()
    for _rackid in RackIDs:
        if _rackid != '-' and _rackid not in Racks:
           Racks[_rackid] = MasterPlate.get(_rackid) 

    # Get Action values
    Actions = dfBC['ACTION'].unique()

    # Create DUMMY Target Rack
    for idx,row in dfBC.iterrows():
        if (row['ACTION'] in  ['NEW']):
            if verbose>0:
                print(f" [{row['ACTION']}] {row['BARCODE']} --> [{row['PLATE_ID']}:{row['WELL_ID']}]")

            # Move from Target to Dummy and add Barcode
            _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
            if _tube_bc.barcode:
                print(f" ERROR [{row['PLATE_ID']}:{row['WELL_ID']}] has existing barcode {_tube_bc.barcode}")                
            _tube_bc.plate_id = DummyRack
            _tube_bc.barcode = row['BARCODE']

            if upload:
                _tube_bc.save()

        elif (row['ACTION'] in ['SAME']):
            if verbose>0:
                print(f" [{row['ACTION']}] --> [{row['PLATE_ID']}:{row['WELL_ID']}]")

            _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
            _tube_bc.plate_id = DummyRack

            if upload:
                _tube_bc.save()

        elif (row['ACTION'] in ['MOVE']):
            if verbose>0:
                print(f" [{row['ACTION']}] {row['BARCODE']} [{row['CURRENT_PLATE_ID']}:{row['CURRENT_WELL_ID']}] --> [{row['PLATE_ID']}:{row['WELL_ID']}]")

            # Move Source Barcodes to DUMMY
            _tube_bc = Racks[row['CURRENT_PLATE_ID']].get_well(row['CURRENT_WELL_ID'])
            _tube_bc.prev_plate_id = str(Racks[row['CURRENT_PLATE_ID']].plate_id)  
            _tube_bc.prev_well_id = row['CURRENT_WELL_ID']  
            _tube_bc.plate_id = DummyRack 
            _tube_bc.well_id = row['WELL_ID']  
            if upload:
                _tube_bc.save()

            # Move Destination Tubes to fill Source Racks, unless same Plate
            if row['PLATE_ID'] != row['CURRENT_PLATE_ID']:
                _tube_mv = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                _tube_mv.plate_id = Racks[row['CURRENT_PLATE_ID']] 
                _tube_mv.well_id = row['CURRENT_WELL_ID']
                if upload:
                    _tube_mv.save()

        elif (row['ACTION'] in ['EMPTY']):
            if verbose>0:
                print(f" [{row['ACTION']}] --> [{row['PLATE_ID']}:{row['WELL_ID']}]")

            _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
            if _tube_bc.barcode:
                print(f" ERROR [{row['PLATE_ID']}:{row['WELL_ID']}] has existing barcode {_tube_bc.barcode}")
            if _tube_bc.n_cmpbatches > 0:
                print(f" ERROR [{row['PLATE_ID']}:{row['WELL_ID']}] has existing compounds {_tube_bc.cmpbatch_id}")               
            _tube_bc.plate_id = DummyRack
            _tube_bc.barcode = None

            if upload:
                _tube_bc.save()

    # Move Dummy Tubes to Target
    DummyRack = MasterPlate.get(DummyRackID)
    if verbose>0:
        print(f" [Save] {TargetRackID}")
    for w in DummyRack.wells:
        DummyRack.wells[w].plate_id = Racks[TargetRackID]
        if upload:
            DummyRack.wells[w].save()



    # # Create New Barcodes
    # if 'NEW' in Actions:
    #     for idx,row in dfBC.iterrows():
    #         if (row['ACTION'] == 'NEW'):
    #             print(f" {row['BARCODE']} ({row['ACTION']})--> [{row['PLATE_ID']}:{row['WELL_ID']}]")
    #             _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
    #             _tube_bc.barcode = row['BARCODE']
    #             if upload:
    #                 _tube_bc.save()

    # # Move existing Barcodes
    # if 'MOVE' in Actions:
    #     # 1st - Move Existing Barcodes to DUMMY and Move Destination Tubes to Source Racks
    #     for idx,row in dfBC.iterrows():
    #         if (row['ACTION'] == 'MOVE'):
    #             print(f" {row['BARCODE']} [{row['CURRENT_PLATE_ID']}:{row['CURRENT_WELL_ID']}] --({row['ACTION']})--> [{row['PLATE_ID']}:{row['WELL_ID']}]")

    #             # Move Source Barcodes to DUMMY
    #             _tube_bc = Racks[row['CURRENT_PLATE_ID']].get_well(row['CURRENT_WELL_ID'])
    #             _tube_bc.prev_plate_id = str(Racks[row['CURRENT_PLATE_ID']].plate_id)  
    #             _tube_bc.prev_well_id = row['CURRENT_WELL_ID']  
    #             _tube_bc.plate_id = DummyRack 
    #             _tube_bc.well_id = row['WELL_ID']  
    #             if upload:
    #                 _tube_bc.save()
                
    #             # Move Destination Tubes to fill Source Racks
    #             _tube_mv = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
    #             _tube_mv.plate_id = Racks[row['CURRENT_PLATE_ID']] 
    #             _tube_mv.well_id = row['CURRENT_WELL_ID']
    #             if upload:
    #                 _tube_mv.save()

    #     # 2nd - Move Tubes from DUMMY to Target Rack           
    #     for idx,row in dfBC.iterrows():
    #         if (row['ACTION'] == 'MOVE'):
    #             # Save with correct PLATE_ID:WELLID
    #             _tube_bc = Racks[row['CURRENT_PLATE_ID']].get_well(row['CURRENT_WELL_ID'])
    #             _tube_bc.plate_id = Racks[row['PLATE_ID']]  
    #             if upload:
    #                 _tube_bc.save()            

        # else:
        #     print(f" {row['BARCODE']} [{row['PLATE_ID']}:{row['WELL_ID']}] --({row['ACTION']})--")

   
    # for idx,row in dfBC.iterrows():
    #     if (row['ACTION'] == 'NEW'):
    #         print(f" {row['BARCODE']} ({row['ACTION']})--> [{row['PLATE_ID']}:{row['WELL_ID']}]")
    #         _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
    #         _tube_bc.barcode = row['BARCODE']
    #         if upload:
    #             _tube_bc.save()

    #     elif (row['ACTION'] != 'SAME'):
    #         print(f" {row['BARCODE']} [{row['CURRENT_PLATE_ID']}:{row['CURRENT_WELL_ID']}] --({row['ACTION']})--> [{row['PLATE_ID']}:{row['WELL_ID']}]")
    #         # Might need to be saved to DummyRack:WELL_ID first - Constraint on PLATE_ID:WELL_ID
    #         _tube_mv = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
    #         _tube_mv.plate_id = DummyRack 
    #         _tube_mv.well_id = row['CURRENT_WELL_ID']
    #         if upload:
    #             _tube_mv.save()
                 
    #         _tube_bc = Racks[row['CURRENT_PLATE_ID']].get_well(row['CURRENT_WELL_ID'])
    #         _tube_bc.prev_plate_id = str(Racks[row['CURRENT_PLATE_ID']].plate_id)  
    #         _tube_bc.prev_well_id = row['CURRENT_WELL_ID']  
    #         _tube_bc.plate_id = Racks[row['PLATE_ID']]  
    #         _tube_bc.well_id = row['WELL_ID']  
    #         if upload:
    #             _tube_bc.save()

    #         # Save with correct PLATE_ID:WELLID
    #         _tube_mv.plate_id = Racks[row['CURRENT_PLATE_ID']]  
    #         if upload:
    #             _tube_mv.save()            

    #     else:
    #         print(f" {row['BARCODE']} [{row['PLATE_ID']}:{row['WELL_ID']}] --({row['ACTION']})--")



