
import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

from dscreen.models import Assay
from dplate.models import Labware, TestPlate, TestWell, MasterPlate, MasterWell
from dorganism.models import Organism, Organism_Batch
from dcell.models import Cell, Cell_Batch

from applib.data.set_fielddata import (set_model_arrayfields, set_model_fields, 
                                       set_model_dicts, set_model_dictarrayfields, 
                                       set_model_fkeys, set_model_from_dict)

#-----------------------------------------------------------------------------
def get_PlatePrep_xlsx(xlsFile, Sheets=[], FillNA='-', **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    PlatePrep_Sheets = {
        'TestPlateList': None,
        'MotherPlates': None,
        'Assays': None,
        'PSPrep': None,
        'HCPrep': None,
        }

    if os.path.isfile(xlsFile):
        fXlsx = open(xlsFile, "rb")
        xls = pd.ExcelFile(fXlsx)
        #xls = pd.ExcelFile(xlsFile)

        if Sheets is None:
            sheets = [PlatePrep_Sheets.keys()]
        for key in Sheets:
            if key in PlatePrep_Sheets:
                PlatePrep_Sheets[key] = pd.read_excel(xls, key)
                PlatePrep_Sheets[key].columns = [c.lower() for c in PlatePrep_Sheets[key].columns]
                if FillNA:
                    PlatePrep_Sheets[key] = PlatePrep_Sheets[key].fillna(FillNA)
            else:
                if valLog:
                    valLog.add_error('Missing Sheet',key,f"XLSX {os.path.basename(xlsFile)}",f"Correct XLSX Sheets {list(PlatePrep_Sheets)}")
        fXlsx.close()
        
    return(PlatePrep_Sheets)

# --------------------------------------------------------------------------------
def read_Motherplates_Prepsheet_XLS(xlFile, prefix=None, **kwargs):
# --------------------------------------------------------------------------------
    PREP_SHEET = 'MotherPlates'

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    _prepSheets = get_PlatePrep_xlsx(xlFile,Sheets=[PREP_SHEET],FillNA=None) 
    # fXlsx = open(xlFile, "rb")
    # xlWB = pd.ExcelFile(fXlsx)

    if _prepSheets[PREP_SHEET] is not None:
        # xDF = xlWB.parse(SheetName)
        # xDF.columns = [c.lower() for c in xDF.columns]

        xDF = _prepSheets[PREP_SHEET]
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

            if valLog:
                if _status == 'New':
                    valLog.add_info("New MotherPlate",
                                    djMP.plate_id, 
                                    f"{djMP.plating}  {djMP.n_wells}w ",
                                    "Select Upload")
                elif _status == 'Exists':
                    valLog.add_warning("MotherPlate Exists ",
                                        djMP.plate_id,
                                        f"{djMP.plating}  {djMP.n_wells}w ",
                                        "Select Overwrite")


            for idx,row in mpwells.iterrows():
                #print(row['compound_id'])
                if not pd.isna(row['compound_id']) and not pd.isna(row['motherwell_id']):
                    # Check if at least compound_id and motherwell_id
                    djWell = djMP.get_well(row['motherwell_id'])
                    set_model_fields(djWell,row,djWell.COPY_FIELDS)
                    set_model_arrayfields(djWell,row,djWell.ARRAY_FIELDS)
                    set_model_dicts(djWell,row,list(djWell.DICTIONARY_FIELDS.keys()),valLog=valLog)
                    set_model_dictarrayfields(djWell,row,djWell.ARRAYDICTIONARY_FIELDS,valLog=valLog)
                    djWell.n_cmpbatches = len(djWell.cmpbatch_lst)
                    
                    validDict = djWell.check_cmpbatch_id()
                    if validDict:
                        validStatus = False
                        dictPl['valid_status'] = False    
                        logger.warning(validDict)
                        row.update(validDict)
                        if valLog:
                            for k in validDict:
                                valLog.add_log(k,'Missing CompoundID',row['compound_id'],f"{djMP.plate_id}:{row['motherwell_id']}",
                                               "Correct PlatePrep or Register Compound")
    
            djMP.add_dilutions(valLog=valLog)
            djMP.set_defaults_model()

            dictPl['new'] = _status == 'New'
            dictPl['plate'] = djMP
            
            lstPl.append(dictPl)
            logger.info(f"[{djMP.plate_id:25s}] - {djMP.plate_type}  {djMP.n_wells}w  [{_status}]")
    else:
        if verbose>0:
            logger.error(f"{PREP_SHEET} not found in {xlFile}")
        if valLog:
            valLog.add_error("Wrong PlatePrep file",
                            f"Sheet: {PREP_SHEET}", 
                            "SheetName not in PlatePrep.XLS",
                            "Correct SheetName")        
    # fXlsx.close()
    return(lstPl)

# --------------------------------------------------------------------------------
def read_TestPlateList_Prepsheet_XLS(xlFile, prefix=None, as_is=False, **kwargs):
# --------------------------------------------------------------------------------
    PREP_SHEETS = ['TestPlateList','Assays']

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    apply_mp = kwargs.get('apply_mp',None)

    logNumbers = {'Processed Assays':0,'New Assays':0, 'Uploaded Assays':0,
                    'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                    'Empty':0}

    _prepSheets = get_PlatePrep_xlsx(xlFile,Sheets=PREP_SHEETS,FillNA='-') 

    # Assays ------------------------------------------------------------------------
    if _prepSheets['Assays'] is not None:

        AssayDict = {}
        ass_Fields = ['assay_subtype','sum_assay_id','test_media', 'test_dye', 'test_enviroment', 'test_time',
                            'test_temperature', 'subculture_type', 'test_addition', ]
        ass_FKeys  = {'organism_id': Organism,'cell_id': Cell}

        validStatus = True

        # Check if correct Sheet
        req_Fields = ['assay_id']
        for _req in req_Fields:
            if _req not in _prepSheets['Assays'].columns:
                validStatus = False
                if valLog:
                    valLog.add_error("Missing Field",
                                    _req, 
                                    "Field not in PlatePrep.XLS [Assays]",
                                    "Correct SheetName")        

        if validStatus:
            logger.info(f"[Assays]")
            for idx,row in _prepSheets['Assays'].iterrows():
                logNumbers['Processed Assays'] += 1
                djAss = Assay.get(row['assay_id'])
                new_assay = False
                validStatus = True
                if djAss is None:
                    djAss = Assay()
                    djAss.assay_id = row['assay_id']
                    if 'organism_id' in row:
                        djAss.assay_type =  row['organism_id']
                    elif 'cell_id' in row:
                        djAss.assay_type =  row['cell_id']    
                    validStatus = set_model_from_dict(djAss,row,list_Fields=ass_Fields, dict_FKeys=ass_FKeys,valLog=valLog)
                    new_assay = True
                    logNumbers['New Assays'] += 1

                # if prgArgs.upload and new_assay and validStatus:
                #     logNumbers['Uploaded Assays'] += 1
                #     djAss.save()

                AssayDict[row['assay_id']] = {'assay': djAss, 'new': new_assay}
            
            valLog.add_info('Assays',f'Listed: {logNumbers["Processed Assays"]}','')
            if logNumbers['New Assays'] > 0:
                valLog.add_warning('Missing Assays',f'New: {logNumbers["New Assays"]}','Assays not registered', 'Update online Assays')

            # logger.info(f"[Assays]: {logNumbers['New Assays']} new assays (of {logNumbers['Processed Assays']}) ")
            # logger.info(f"[Assays]: New assays uploaded {logNumbers['Uploaded Assays']} [Upload: {prgArgs.upload}]")
            # logger.info(f"[Assays]")


    # TestPlates ------------------------------------------------------------------------
    if _prepSheets['TestPlateList'] is not None:

        TestPlateDict = {}
        tp_Fields = ['plating','test_media', 'test_dye','test_additive', 'processing', 'issues','control_layout' ]
        tp_FKeys  = {'assay_id': Assay,'test_cellbatch_id': Cell_Batch, 'test_orgbatch_id': Organism_Batch, 'labware_id': Labware}
        tp_Dicts = ['result_type']
        tp_Arrays = {'motherplate_ids':['motherplate_id','motherplate2_id'],'synergy_cmpbatches':['syn_compounds_ab','syn_compounds_pot']}
        tp_FKeyArrays = {'motherplate_ids':{'model':MasterPlate, 'fields': ['motherplate_id','motherplate2_id']},
                            'synergy_cmpbatches':{'model':Organism_Batch, 'fields':['syn_compounds_ab','syn_compounds_pot']}
                        }
        validStatus = True

        # Check if correct Sheet
        req_Fields = ['testplate_id','result_type','assay','assay_id','test_strain','test_cell','control_layout','motherplate_id']
        for _req in req_Fields:
            if _req not in _prepSheets['TestPlateList'].columns:
                validStatus = False
                if valLog:
                    valLog.add_error("Missing Field",
                                    _req, 
                                    "Field not in PlatePrep.XLS [TestPlateList]",
                                    "Correct SheetName")        

        if validStatus:
            logger.info(f"[TestPlateList]")
            _prepSheets['TestPlateList'].rename(columns={"test_strain": "test_orgbatch_id", 
                                                        "test_cell": "test_cellbatch_id"},inplace=True)
            for idx,row in _prepSheets['TestPlateList'].iterrows():
                logNumbers['Processed Plates'] += 1
                validStatus = True
                djTP = TestPlate.get(row['testplate_id'],WellData=apply_mp)
                if djTP is None:
                    valLog.add_error(f"MissingTestPlate", row['testplate_id'],'Testplate not found','Upload ReadOuts or Check TestPlate_ID')
                    logNumbers['New Plates'] += 1
                else:
                    validStatus = set_model_from_dict(djTP,row,
                                                      list_Fields=tp_Fields, 
                                                      dict_Arrays=tp_Arrays, 
                                                      list_Dicts=tp_Dicts, 
                                                      dict_FKeys=tp_FKeys,
                                                      dict_FKeyArrays= tp_FKeyArrays,
                                                      valLog=valLog)
                    
                    TestPlateDict[row['testplate_id']] = {'plate':djTP}

            valLog.add_info('TestPlates',f'Listed: {logNumbers["Processed Plates"]}','Testplate listed',)

    return TestPlateDict,AssayDict