import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime
from tqdm import tqdm

import logging
logger = logging.getLogger(__name__)

from apputil.models import Dictionary
from dplate.models import MasterPlate, MasterWell, TestPlate, TestWell
from dsample.models import Compound_Batch,COADD_Compound
#from applib.data.set_fielddata import set_model_arrayfields, set_model_fields, set_model_dicts, set_model_fkeys, set_model_dictarrayfields
from decimal import Decimal


#--------------------------------------------------------------------------------------------------------------
# Add MotherWell to TestWell
#--------------------------------------------------------------------------------------------------------------
def add_mother_to_testwell(djMW,djTW):
    
    FIELD_LIST = [['cmpbatch_lst','cmpbatch_lst'],
                  ['set_lst','set_lst'],
                  ['test_conc_lst','conc_lst'],
                  ['test_conc_unit_lst','conc_unit_lst'],
                  ['test_conc_type_lst','conc_type_lst'],
                ]
        
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
def add_mother_to_testplate(djMP,djTP):
    if hasattr(djTP,'wells') and hasattr(djMP,'wells'):
        if djTP.n_wells == djMP.n_wells:
            for w in djMP.wells:
                add_mother_to_testwell(djMP.wells[w],djTP.wells[w])
        else:
            logger.info(f" [MP->TP] Different n_wells {djMP.n_wells} -> {djTP.n_wells}")     
    else:
        logger.info(f" [MP->TP] Missing well data {djMP.plate_id} -> {djTP.plate_id}")
        
        
#--------------------------------------------------------------------------------------------------------------
# Read CpOz Plate Files
#--------------------------------------------------------------------------------------------------------------
def read_cpoz_edr(EDR_File, Skip_File=None, PlateSize=384, TestVolume=50, TestConc='uM'): 
    
    # 3P_UQA37_CKL01_OASCAF_8ptEDR_C3-4_EDR Export
    # QCL_Sample_Number	Destination_Plate_Barcode	Destination_Row	Destination_Column	RowCol	Destination_Concentration	Actual_Volume	Backfill_Volume	Total_Volume

    SN_IDs = {}
    TPltDict = {}
    SNidDict = {}
    nSmpDict = {}
     
    if os.path.exists(EDR_File):
        EDR_Base = os.path.basename(EDR_File)
        logger.info(f" [CpOz EDR  ] {EDR_File}")
        pEDR = pd.read_csv(EDR_File,sep='\t')
                
        TPltLst = pEDR['Destination_Plate_Barcode'].unique()
        SNidLst = pEDR['QCL_Sample_Number'].unique()

        # Get Testplate instances --------------------------
        for PID in TPltLst:
            TPltDict[PID] = TestPlate.get(PID)
            if TPltDict[PID] is None:
                logger.info(f" [CpOz EDR  ] TestPlate {PID} [New]")
                TPltDict[PID] = TestPlate.new(PID,PlateSize,WellData=True)            
            else:
                logger.info(f" [CpOz EDR  ] TestPlate {PID} [Exists]")
            nSmpDict[PID] = 0

        # Get Cmpound instances ----------------------------
        ValidStatus = True
        for SID in tqdm(SNidLst,desc='QCL SN'):
            if 'DMSO' not in SID:
                if COADD_Compound.objects.filter(cpoz_sn=SID).exists():
                    qryCmp = COADD_Compound.objects.filter(cpoz_sn=SID)
                    SNidDict[SID] = qryCmp[0].cmpbatch_id
                    #logger.info(f" [CpOz EDR  ] QCL Sample {SID} : { SNidDict[SID]}")
                else:
                    logger.error(f" [CpOz EDR  ] QCL Sample {SID} Not Found")
                    ValidStatus = False
        
        # Get TestWells ----------------------------
        if ValidStatus:
            djPCT = Dictionary.get(TestWell.DICTIONARY_FIELDS['solvent_conc_unit'],'pct')
            for idx,row in tqdm(pEDR.iterrows(), total=len(pEDR), desc='Wells'):
                SID = row['QCL_Sample_Number']
                if 'DMSO' not in SID:
                    PID = row['Destination_Plate_Barcode']
                    WID = row['RowCol']
                    TPltDict[PID].wells[WID].cmpbatch_id = SNidDict[SID]
                    TPltDict[PID].wells[WID].cmpbatch_lst = [SNidDict[SID]]
                    TPltDict[PID].wells[WID].n_cmpbatches =1
                                     
                    tConc = row['Destination_Concentration'] * 1000 * 1000
                    sConc = 100* (row['Total_Volume'])/(TestVolume *1000)
                    
                    TPltDict[PID].wells[WID].conc_lst = [tConc]
                    TPltDict[PID].wells[WID].conc_unit_lst = [TestConc]
                    TPltDict[PID].wells[WID].solvent = 'DMSO'
                    TPltDict[PID].wells[WID].solvent_conc = sConc
                    TPltDict[PID].wells[WID].solvent_conc_unit = djPCT
                    nSmpDict[PID] += 1
        
        # Add Plate Infos
        for PID in TPltLst:
            TPltDict[PID].n_samples = nSmpDict[PID]

        # Read Skip File
        if Skip_File is not None:
            if os.path.exists(Skip_File):
                logger.info(f" [CpOz Skips] {Skip_File}")
                
                pSkips = pd.read_csv(Skip_File,sep='\t')
                for idx,row in pSkips.iterrows():
                    PID = row['Destination_Plate_Barcode']
                    if PID in TPltDict:
                        WID = None
                        if 'Tgt Coord' in row:
                            WID = row['Tgt Coord']
                        elif 'Destination_Row' in row:
                            WID = TPltDict[PID].well_pos((row['Destination_Row'],row['Destination_Column']))
                        if WID:
                            TPltDict[PID].wells[WID].is_skip = True
                            logger.info(f" [CpOz Skips] {PID} {WID}")

    return(TPltDict)