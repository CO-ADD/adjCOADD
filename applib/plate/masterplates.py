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
class Barcode_Storage():

    def __init__(self,**kwargs):

        self.rack_size = 96
        self.rack_type = 'Storage'
        self.rack_labware = 'FLUIDX_10mL'
        self.storage_size = 400
        self.storage_name = "ORPHAN_BARCODES"
        self.barcode_location = {}
        self.dummy_storage_name = "DUMMY_BARCODES"
        self.dummy_rack_name = "DUMMY_RACK"

        self.verbose = kwargs.get('verbose',0)

        self.dummy_storage = self._load_or_create_dummy_rack(self.dummy_storage_name, self.storage_size, reset=True)
        self.dummy_storage_id = str(self.dummy_storage.plate_id)
        self.dummy_rack = self._load_or_create_dummy_rack(self.dummy_rack_name, self.rack_size, reset=True)
        self.dummy_rack_id = str(self.dummy_rack.plate_id)

        self._load_or_create_barcode_storage()

        if self.verbose > 0:
            print(f" [Barcode Storage] {self.storage_id} : {len(self.barcode_location)} barcodes")
            print(f" [Dummy   Storage] {self.dummy_storage_id} ")
            print(f" [Dummy   Rack   ] {self.dummy_rack_id} : ")


    # -----------------------------------------
    def _load_or_create_barcode_storage(self):
    # -----------------------------------------
        self.storage_id = f"{self.storage_name}_{self.storage_size}"
        if MasterPlate.exists(self.storage_id):
            self.storage = MasterPlate.get(self.storage_id)
            self._get_stored_barcode_locations()
        else:
            self.storage = MasterPlate().new(self.storage_id,self.storage_size,self.rack_type, WellData=False)
            self.storage.set_defaults_model()
            self.storage.save()
            self.barcode_location = {}

    # -----------------------------------------
    def _load_or_create_dummy_rack(self,rack_name, rack_size, reset=True):
    # -----------------------------------------
        _rack_id = f"{rack_name}_{rack_size}"
        if MasterPlate.exists(_rack_id):
            _dummy= MasterPlate.get(_rack_id, FillMissing=False)
            if reset:
                _dummy.delete_wells() 
        else:
            _dummy = MasterPlate().new(_rack_id,rack_size,self.rack_type, WellData=False)
            _dummy.set_defaults_model()
            _dummy.save()
            # Set aStatus -9 : Invisible
            _dummy.remove()
        return(_dummy)

    # # -----------------------------------------
    # def _load_or_create_dummy_storage(self):
    # # -----------------------------------------
    #     self.dummy_storage_id = f"{self.dummy_name}_{self.storage_size}"
    #     if MasterPlate.exists(self.dummy_storage_id):
    #         self.dummy_storage= MasterPlate.get(self.dummy_storage_id)
    #     else:
    #         self.dummy_storage = MasterPlate().new(self.dummy_storage_id,self.storage_size,self.rack_type, WellData=False)
    #         self.dummy_storage.set_defaults_model()
    #         self.dummy_storage.save()

    # -----------------------------------------
    def _get_stored_barcode_locations(self):
    # -----------------------------------------
        self.barcode_location = {}
        for w in self.storage.wells:
            if self.storage.wells[w].barcode:
                self.barcode_location[w] = self.storage.wells[w].barcode

    # -----------------------------------------
    def _next_location(self):
    # -----------------------------------------
        for w in self.storage.wells:
            if w not in self.barcode_location:
                return(w)
        return(None)

    # -----------------------------------------
    def _reload_barcode_storage_dummys(self):
    # -----------------------------------------
        self._load_or_create_barcode_storage()
        self.dummy_storage = self._load_or_create_dummy_rack(self.dummy_storage_name, self.storage_size, reset=True)
        self.dummy_rack = self._load_or_create_dummy_rack(self.dummy_rack_name, self.rack_size, reset=True)

    # -----------------------------------------
    def store_tube(self,Tube, RemoveEmpty=False):
        if Tube.barcode:
            _nw = self._next_location()

            _t_plate_id = Tube.plate_id
            _t_well_id  = Tube.well_id

            # Move empty STORE Tube to DUMMY Tube
            #self.storage.wells[_nw].plate_id = self.dummy_storage
            #self.storage.wells[_nw].save()

            # Move SOURCE Tube to STORE 
            Tube.prev_plate_id = str(_t_plate_id)
            Tube.prev_well_id = _t_well_id
            Tube.plate_id = self.storage
            Tube.well_id = _nw
            Tube.save()
            self.barcode_location[_nw] = Tube.barcode

        elif RemoveEmpty:
            Tube.delete()

            # Move empty DUMMY Tube to SOURCE
            #self.storage.wells[_nw].plate_id = _t_plate_id
            #self.storage.wells[_nw].well_id  = _t_well_id
            #self.storage.wells[_nw].save()

    # -----------------------------------------
    def store_rack(self,Rack):
        for w in Rack.wells:
            self.store_tube(Rack.wells[w])

    # -----------------------------------------
    def move_tube(self,Tube, PlateID, WellID, RemoveEmpty=False):
        if Tube.barcode:
            # Move Barcode to new location 
            Tube.plate_id = PlateID
            Tube.well_id = WellID
            Tube.save()

        elif RemoveEmpty:
            Tube.delete()
    # -----------------------------------------
    def _read_barcode_scan(self,csvFile):    
        self.barcodes = pd.read_csv(csvFile)
        self.barcodes.columns = ['PLATE_ID', 'WELL_ID', 'BARCODE']


    #--------------------------------
    @staticmethod
    def _apply_current_barcode_location(s):
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
                elif (s['CURRENT_PLATE_ID'] == s['PLATE_ID']):
                    s['ACTION'] = 'SWAP'
                else:
                    s['ACTION'] = 'MOVE'
            else:
                s['CURRENT_PLATE_ID'] = '-'
                s['CURRENT_WELL_ID'] = '-'
                s['ACTION'] = 'NEW'
        return(s)

    # -------------------------------------------------------------------------
    def update_barcode_location(self,csvFile, upload=False):
    # -------------------------------------------------------------------------

        self._read_barcode_scan(csvFile)

        # Get/Create Target Rack ----------------------------------------
        targetRackID = self.barcodes['PLATE_ID'].unique()[0]
        if MasterPlate.exists(targetRackID):
            # Move existing Barcodes into Storage and make targetRack empty 
            _targetRack = MasterPlate.get(targetRackID)
            self.store_rack(_targetRack)

            Racks = {targetRackID:MasterPlate.get(targetRackID)}            
        else:
            Racks = {targetRackID:MasterPlate().new(targetRackID,self.rack_size,self.rack_type, WellData=True)}
            Racks[targetRackID].labware_id = Labware.get(self.rack_labware)
            Racks[targetRackID].save()


        print("-------------------------")
        # Get all Source Racks ----------------------------------------
        self.barcodes = self.barcodes.apply(self._apply_current_barcode_location,axis=1)
        RackIDs = self.barcodes['CURRENT_PLATE_ID'].unique()
        for _rackid in RackIDs:
            if _rackid != '-' and _rackid not in Racks:
                Racks[_rackid] = MasterPlate.get(_rackid)
                print(f" Load RackID {_rackid}") 

        # Get Action values
        Actions = self.barcodes['ACTION'].unique()

        print("-------------------------")
        # Create TargetRack in DUMMY
        for idx,row in self.barcodes.iterrows():

            if (row['ACTION'] in  ['MOVE','SWAP','SAME']):
 
                # Move CURRENT Barcodes to DUMMY Rack -----------------------------------
                _tube_bc = Racks[row['CURRENT_PLATE_ID']].get_well(row['CURRENT_WELL_ID'])
                _tube_bc.plate_id = self.dummy_rack 
                _tube_bc.well_id = row['WELL_ID']  

                # Set SOURCE information into prev_plate_id/well_id
                #     STORE tube should have prev_plate_id/well_id
                if row['CURRENT_PLATE_ID'] != self.storage_id:
                    _tube_bc.prev_plate_id = str(Racks[row['CURRENT_PLATE_ID']].plate_id)  
                    _tube_bc.prev_well_id = row['CURRENT_WELL_ID']  
                else:
                    self.barcode_location.pop(row['CURRENT_WELL_ID'],None)

                _tube_bc.save()
                self.dummy_rack.wells[_tube_bc.well_id] = _tube_bc
                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} [{row['CURRENT_PLATE_ID']}:{row['CURRENT_WELL_ID']}] --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")

                # Move non-empty Target Tube to CURRENT Rack ---------------------------------
                _tube_mv = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])

                # self.move_tube(_tube_mv,dest_plate,dest_well,RemoveEmpty=True)
                #
                #
                if _tube_mv.barcode:
                    _tube_mv.plate_id = Racks[row['CURRENT_PLATE_ID']] 
                    _tube_mv.well_id = row['CURRENT_WELL_ID']                
                    _tube_mv.save()
                    print(f" TARGET->CURRENT   {row['PLATE_ID']} {row['WELL_ID']} -> [{_tube_mv.id} {str(_tube_mv.plate_id)} {_tube_mv.well_id} {_tube_mv.barcode}]")
                else:
                    # Remove empty Target Tube
                    _tube_mv.delete()

            elif (row['ACTION'] in  ['NEW']):

                 # Move CURRENT Barcodes to DUMMY Rack -----------------------------------
                _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                _tube_bc.plate_id = self.dummy_rack
                _tube_bc.barcode = row['BARCODE']
                _tube_bc.save()
                self.dummy_rack.wells[_tube_bc.well_id] = _tube_bc
                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")

            elif (row['ACTION'] in  ['EMPTY']):

                # Option: 
                # Remove from Racks[row['PLATE_ID']] if empty, otherwise move to STORAGE
                # Leave DUMMY empty
                # 
                # self.store_tube(tube_bc,RemoveEmpty=True)


                # Move CURRENT Barcodes to DUMMY Rack -----------------------------------
                _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                _tube_bc.plate_id = self.dummy_rack
                _tube_bc.save()
                self.dummy_rack.wells[_tube_bc.well_id] = _tube_bc
                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")

        print("-------------------------")
        for w in self.dummy_rack.wells:
            print(f" {w} [{self.dummy_rack.wells[w]}]")

        print("-------------------------")
        # Move DUMMY to Target Rack
        for idx,row in self.barcodes.iterrows():
            # Move DUMMY Barcodes to TARGET Rack

            _tube_bc=self.dummy_rack.get_well(row['WELL_ID'])
            _px = str(_tube_bc.plate_id)
            _wx = str(_tube_bc.well_id)

            _tube_bc.plate_id = Racks[row['PLATE_ID']] 
            _tube_bc.well_id = row['WELL_ID']  
            _tube_bc.save()

            print(f" DUMMY->TARGET   {_px} {_wx} -> {str(_tube_bc.plate_id)} {_tube_bc.well_id} {_tube_bc.barcode}")

        # Reset Barcode Storage
        self._reload_barcode_storage_dummys()

# # --------------------------------------------------------------------------------
# def read_barcode_csv(csvFile, add_current_location=True, **kwargs):
# # --------------------------------------------------------------------------------

#     def apply_current_location(s):
#         if s['BARCODE'] == 'NO READ':
#             s['CURRENT_PLATE_ID'] = '-'
#             s['CURRENT_WELL_ID'] = '-'
#             s['ACTION'] = 'EMPTY'
#         else:
#             _tube = MasterWell.get(None,None,s['BARCODE'])
#             if _tube:
#                 s['CURRENT_PLATE_ID'] = str(_tube.plate_id)
#                 s['CURRENT_WELL_ID'] = str(_tube.well_id)

#                 if ((s['CURRENT_PLATE_ID'] == s['PLATE_ID']) and (s['CURRENT_WELL_ID'] == s['WELL_ID'])):
#                     s['ACTION'] = 'SAME'
#                 elif (s['CURRENT_PLATE_ID'] == s['PLATE_ID']):
#                     s['ACTION'] = 'SWAP'
#                 else:
#                     s['ACTION'] = 'MOVE'
#             else:
#                 s['CURRENT_PLATE_ID'] = '-'
#                 s['CURRENT_WELL_ID'] = '-'
#                 s['ACTION'] = 'NEW'
#         return(s)

#     dfBarcode = pd.read_csv(csvFile)
#     dfBarcode.columns = ['PLATE_ID', 'WELL_ID', 'BARCODE']
#     if add_current_location:
#         dfBarcode = dfBarcode.apply(apply_current_location,axis=1)
    
#     return(dfBarcode)             


# # --------------------------------------------------------------------------------
# def get_dummy_rack(DummyID='MP_DUMMY_MOVE', PlateSize=96, LabewareID='FLUIDX_10mL',upload=False):
# # --------------------------------------------------------------------------------
#     RACK_TYPE = 'Storage'
#     _plate_id = f"{DummyID}_{PlateSize}"
#     if MasterPlate.exists(_plate_id):
#         return(MasterPlate.get(_plate_id))
#     else:
#         _MP = MasterPlate().new(_plate_id,PlateSize,RACK_TYPE, WellData=False)
#         _MP.labware_id = Labware.get(LabewareID)
#         if upload:
#             _MP.set_defaults_model()
#             _MP.save()
#             _MP.remove()
#         return(_MP)

# # --------------------------------------------------------------------------------
# def move_rack(RackID,DummyID='MP_DUMMY_COPY', PlateSize=96, LabewareID='FLUIDX_10mL',upload=False):
# # --------------------------------------------------------------------------------
#     RACK_TYPE = 'Storage'
#     _plate_id = f"{DummyID}_{PlateSize}"
    

# --------------------------------------------------------------------------------
def update_barcode_location(csvFile, LabewareID = 'FLUIDX_10mL', upload=False, verbose=0, debug_step=0):
# --------------------------------------------------------------------------------
    RACK_SIZE = 96
    RACK_TYPE = 'Storage'

    #--------------------------------
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
                elif (s['CURRENT_PLATE_ID'] == s['PLATE_ID']):
                    s['ACTION'] = 'SWAP'
                else:
                    s['ACTION'] = 'MOVE'
            else:
                s['CURRENT_PLATE_ID'] = '-'
                s['CURRENT_WELL_ID'] = '-'
                s['ACTION'] = 'NEW'
        return(s)
    #--------------------------------

    # Get/Create Barcode Storage
    bcStorage = Barcode_Storage()

    # Setup DummyRack
    dummyRack = get_dummy_rack(upload=upload)
    dummyRackID = str(dummyRack.plate_id)

    # Read Barcode Scan
    dfBC = pd.read_csv(csvFile)
    dfBC.columns = ['PLATE_ID', 'WELL_ID', 'BARCODE']

    # Get and store, or Create new Target Rack --------------------------------------
    targetRackID = dfBC['PLATE_ID'].unique()[0]
    if MasterPlate.exists(targetRackID):
        Racks = {targetRackID:MasterPlate.get(targetRackID)}
        target_is_empty = False

        # Move existing Barcodes into Storage and make targetRack empty
        bcStorage.store_rack(Racks[targetRackID])
    else:
        Racks = {targetRackID:MasterPlate().new(targetRackID,RACK_SIZE,RACK_TYPE, WellData=True)}
        Racks[targetRackID].labware_id = Labware.get(LabewareID)
        target_is_empty = True
        if upload:
            Racks[targetRackID].save()

    # Get all Source Racks ----------------------------------------
    RackIDs = dfBC['CURRENT_PLATE_ID'].unique()
    for _rackid in RackIDs:
        if _rackid != '-' and _rackid not in Racks:
           Racks[_rackid] = MasterPlate.get(_rackid) 

    # Get Action values
    dfBC = dfBC.apply(apply_current_location,axis=1)
    Actions = dfBC['ACTION'].unique()

    # Move Existing Barcodes to DummyRack
    for idx,row in dfBC.iterrows():
        if (row['ACTION'] in  ['SAME','SWAP']):
            if verbose>0:
                print(f" [{row['ACTION']}] --> [{row['PLATE_ID']}:{row['WELL_ID']}]")

            _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
            _tube_bc.plate_id = DummyRack

            if upload:
                _tube_bc.save()

        if (row['ACTION'] in  ['MOVE']):
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
        
# --------------------------------------------------------------------------------
def x_update_barcode_location(csvFile, LabewareID = 'FLUIDX_10mL', upload=False, verbose=0, debug_step=0):
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
    if debug_step > 1:  
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



