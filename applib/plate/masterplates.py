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
    """
    Barcode Handling Class
        Racks (MasterPlates as 'Storage') and Tubes (MasterWell with Barcodes) 
        Use Barcode Scan to Assign new Barcodes, to Move existing Barcodes
        
        Uses a ORPHAN_BARCODES rack for Tubes without a Rack Location
            Might need cleaning if many Barcodes are physically removed, or size increased

    """
    def __init__(self, **kwargs):

        self.rack_size = 96
        self.rack_type = 'Storage'
        self.rack_labware = 'FLUIDX_10mL'

        self.storage_size = 2000
        self.storage_id = "ORPHAN_BARCODES"
        self.storage_barcodes = {}

        self.dummy_storage_name = "DUMMY_BARCODES"
        self.dummy_rack_name = "DUMMY_RACK"

        # ----------------------------------------
        self.verbose = kwargs.get('verbose',0)

        self.dummy_rack = self._load_or_create_dummy_rack(self.dummy_rack_name, self.rack_size, reset=True)
        self.dummy_rack_id = str(self.dummy_rack.plate_id)

        self._load_or_create_barcode_storage()

        if self.verbose > 0:
            print(f" [Barcode Storage] {self.storage_id} : {len(self.storage_barcodes)} barcodes")
            print(f" [Dummy   Rack   ] {self.dummy_rack_id} : ")


    # -----------------------------------------
    def _load_or_create_barcode_storage(self):
    # -----------------------------------------
        if MasterPlate.exists(self.storage_id):
            self.storage = MasterPlate.get(self.storage_id)
            self._get_stored_barcode_locations()
        else:
            self.storage = MasterPlate().new(self.storage_id,self.storage_size,self.rack_type, WellData=False)
            self.storage.set_defaults_model()
            self.storage.save()
            self.storage_barcodes = {}

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

    # -----------------------------------------
    def _get_stored_barcode_locations(self):
    # -----------------------------------------
        self.storage_barcodes = {}
        for w in self.storage.wells:
            if self.storage.wells[w].barcode:
                self.storage_barcodes[w] = self.storage.wells[w].barcode

    # -----------------------------------------
    def _next_location(self):
    # -----------------------------------------
        for w in self.storage.wells:
            if w not in self.storage_barcodes:
                return(w)
        return(None)

    # -----------------------------------------
    def _reload_barcode_storage_dummys(self):
    # -----------------------------------------
        self._load_or_create_barcode_storage()
        self.dummy_rack = self._load_or_create_dummy_rack(self.dummy_rack_name, self.rack_size, reset=True)

    # -----------------------------------------
    def store_tube(self,Tube, RemoveEmpty=False):
        if Tube:
            if Tube.barcode:
                _nw = self._next_location()

                # Move SOURCE Tube to STORE 
                Tube.prev_plate_id = str(Tube.plate_id)
                Tube.prev_well_id = Tube.well_id
                Tube.plate_id = self.storage
                Tube.well_id = _nw
                Tube.save()
                self.storage_barcodes[_nw] = Tube.barcode

            elif RemoveEmpty:
                if Tube.id:
                    Tube.delete()

    # -----------------------------------------
    def store_rack(self,Rack):
        for w in Rack.wells:
            self.store_tube(Rack.wells[w])

    # -----------------------------------------
    def move_tube(self,Tube, PlateID, WellID, RemoveEmpty=False):
        if Tube:
            if Tube.barcode:
                # Move Barcode to new location 
                Tube.plate_id = PlateID
                Tube.well_id = WellID
                Tube.save()

            elif RemoveEmpty:
                if Tube.id:
                    Tube.delete()
    # -----------------------------------------
    def _read_barcode_scan(self,csvFile):    
        self.barcodes = pd.read_csv(csvFile)
        self.barcodes.columns = ['PLATE_ID', 'WELL_ID', 'BARCODE']

    #--------------------------------
    @staticmethod
    def _apply_current_barcode_location(s):
        if s['BARCODE'] == 'NO READ':
            s['SOURCE_PLATE_ID'] = '-'
            s['SOURCE_WELL_ID'] = '-'
            s['ACTION'] = 'EMPTY'
        else:
            _tube = MasterWell.get(None,None,s['BARCODE'])
            if _tube:
                s['SOURCE_PLATE_ID'] = str(_tube.plate_id)
                s['SOURCE_WELL_ID'] = str(_tube.well_id)

                if ((s['SOURCE_PLATE_ID'] == s['PLATE_ID']) and (s['SOURCE_WELL_ID'] == s['WELL_ID'])):
                    s['ACTION'] = 'SAME'
                elif (s['SOURCE_PLATE_ID'] == s['PLATE_ID']):
                    s['ACTION'] = 'SWAP'
                else:
                    s['ACTION'] = 'MOVE'
            else:
                s['SOURCE_PLATE_ID'] = '-'
                s['SOURCE_WELL_ID'] = '-'
                s['ACTION'] = 'NEW'
        return(s)

    # -------------------------------------------------------------------------
    def update_barcode_location(self, csvFile, upload=False):
    # -------------------------------------------------------------------------

        self._read_barcode_scan(csvFile)
        # Get/Create TARGET Rack ----------------------------------------
        targetRackID = self.barcodes['PLATE_ID'].unique()[0]
        print(f" [Update Barcode Location] {targetRackID} : {csvFile} ")

        if MasterPlate.exists(targetRackID):
            # Move existing Barcodes into STORE 
            _targetRack = MasterPlate.get(targetRackID)
            self.store_rack(_targetRack)

            # Reload the empty TARGET Rack
            Racks = {targetRackID:MasterPlate.get(targetRackID)}            
        else:
            # Create the empty TARGET Rack
            Racks = {targetRackID:MasterPlate().new(targetRackID,self.rack_size,self.rack_type, WellData=True)}
            Racks[targetRackID].labware_id = Labware.get(self.rack_labware)
            Racks[targetRackID].save()


        # Get all SOURCE Racks ----------------------------------------
        self.barcodes = self.barcodes.apply(self._apply_current_barcode_location,axis=1)
        RackIDs = self.barcodes['SOURCE_PLATE_ID'].unique()
        for _rackid in RackIDs:
            if _rackid != '-' and _rackid not in Racks:
                Racks[_rackid] = MasterPlate.get(_rackid)

        # Create DUMMY as Target
        for idx,row in self.barcodes.iterrows():

            if (row['ACTION'] in  ['MOVE','SWAP','SAME']): 
                # Move SOURCE Barcodes to DUMMY Rack -----------------------------------
                _tube_bc = Racks[row['SOURCE_PLATE_ID']].get_well(row['SOURCE_WELL_ID'])
                _tube_bc.plate_id = self.dummy_rack 
                _tube_bc.well_id = row['WELL_ID']  

                # Set SOURCE information into prev_plate_id/well_id
                if row['SOURCE_PLATE_ID'] != self.storage_id:
                    _tube_bc.prev_plate_id = str(Racks[row['SOURCE_PLATE_ID']].plate_id)  
                    _tube_bc.prev_well_id = row['SOURCE_WELL_ID']  
                else:
                    # STORE tube should have prev_plate_id/well_id already
                    self.storage_barcodes.pop(row['SOURCE_WELL_ID'],None)

                _tube_bc.save()
                self.dummy_rack.wells[_tube_bc.well_id] = _tube_bc

                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} [{row['SOURCE_PLATE_ID']}:{row['SOURCE_WELL_ID']}] --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")

                # Move non-empty Target Tube to SOURCE Rack, remove empty Tube ----------
                _tube_mv = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                self.move_tube(_tube_mv,Racks[row['SOURCE_PLATE_ID']],row['SOURCE_WELL_ID'],RemoveEmpty=True)

            elif (row['ACTION'] in  ['NEW']):
                # Add Barcode to TARGET and move to DUMMY Rack -----------------------------------
                _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                _tube_bc.plate_id = self.dummy_rack
                _tube_bc.barcode = row['BARCODE']
                _tube_bc.save()
                self.dummy_rack.wells[_tube_bc.well_id] = _tube_bc
                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")

            elif (row['ACTION'] in  ['EMPTY']):
                # Remove from Racks[row['PLATE_ID']] if empty, otherwise move to STORAGE
                _tube_bc = Racks[row['PLATE_ID']].get_well(row['WELL_ID'])
                if self.verbose>0:
                    print(f" [{row['ACTION']}] {row['BARCODE']} --> [{_tube_bc.id} {str(_tube_bc.plate_id)}:{_tube_bc.well_id} {_tube_bc.barcode}]")
                self.store_tube(_tube_bc,RemoveEmpty=True)

        # if self.verbose>0:
        #     for w in self.dummy_rack.wells:
        #         print(f" {w} [{self.dummy_rack.wells[w]}]")

        # Move DUMMY to Target Rack
        for idx,row in self.barcodes.iterrows():
            # Move DUMMY Barcodes to TARGET Rack
            _tube_bc=self.dummy_rack.get_well(row['WELL_ID'])
            self.move_tube(_tube_bc,Racks[row['PLATE_ID']],row['WELL_ID'],RemoveEmpty=True)

        # Reset Barcode Storage
        self._reload_barcode_storage_dummys()
        if len(self.storage_barcodes) > 0:
            print(f" [Barcode Storage] {self.storage_id} : {len(self.storage_barcodes)} barcodes")


