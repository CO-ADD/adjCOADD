
import os,sys
import csv
import re
import numpy as np
import pandas as pd
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)

from dplate.models import TestPlate
from decimal import Decimal

# --------------------------------------------------------------------------------
def multimodereader_xls(xlFile, prefix=None, as_is=False):
# --------------------------------------------------------------------------------

    xlWB = pd.ExcelFile(xlFile)
    lstPl = []
    for xSheet in xlWB.sheet_names:
        xDF = xlWB.parse(xSheet, header = None)

        if len(xDF)>0:
            djTP = None
            if xDF[0][0] == "Application: Tecan i-control":
                djTP = read_iControl_xlsheet(xSheet,xDF,prefix=prefix)
                
            # elif xDF[0][0] == "Experiment" or xDF[0][1] == "Experiment":
            #     xPl = readPlate_Gen5_sheet(xSheet,xDF,prefix=prefix)
            
            # elif "CLARIOstar" in xDF[0][3]:
            #     xPl = readPlate_BMG_sheet(xSheet,xDF,prefix=prefix)

            if djTP:
                djTP.input_file = os.path.split(xlFile)[1]

                lstPl.append(djTP)
                logger.info(f"[{djTP.plate_id:25s}] - {djTP.reader}  {djTP.n_wells} {djTP.readout_type}")
            else:
                logger.info(f"[{xSheet}] - Unknown PlateReader Format")
    return(lstPl)


#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - - Tecan iControl XLSX.SHEET
#--------------------------------------------------------------------------------------------------------------
def read_iControl_xlsheet(xSheet,xDF,prefix=None):
    max_PlateID_len = 25

    plateDict = {}
    plateDict['READS'] = []

    paramLst = ['Mode','System','User','Plate-ID','Barcode',
                'Measurement Wavelength',
                'Excitation Wavelength','Emission Wavelength','Excitation Bandwidth','Emission Bandwidth',
                'Gain','Number of Flashes','Flash Frequency','Integration Time','Lag Time',
                'Settle Time'
                ]
    nWells = 0
    _flMatrix = 0
    _matrix = []
    _ncol = 0

    for r in range(len(xDF)):
        #print(xDF[r])
        for k in paramLst:
            if k in str(xDF[0][r]):
                if str(xDF[4][r]) != "nan":
                    if k.upper() in plateDict:
                        plateDict[k.upper()] += ";"+str(xDF[4][r])
                    else:
                        plateDict[k.upper()] = str(xDF[4][r])

        if 'End Time:' in str(xDF[0][r]):
            try:
                plateDict['TEST_DATE'] = datetime.strptime(str(xDF[1][r]), '%d/%m/%Y %H:%M:%S %p')
            except:
                plateDict['TEST_DATE'] = datetime.strptime(str(xDF[1][r]), '%m/%d/%Y %H:%M:%S %p')

        # MatrixRead - Row 
        if _flMatrix > 0 :
            _row = []
            if str(xDF[0][r]) in TestPlate.ROW_LABELS:
                for c in range(_ncol):
                    _row.append(xDF[c+1][r])
                _matrix.append(_row)
            else:
                plateDict['READS'].append(_matrix)
                _flMatrix = 0
                _matrix = []

        # MatrixRead - Init 
        if '<>' in str(xDF[0][r]):
            if xDF[12][r] == 12:
                nWells = 96
                _ncol = 12
                _nrow = 8
                _flMatrix = 1
            if xDF[24][r] == 24:
                nWells = 384
                _ncol = 24
                _nrow = 16
                _flMatrix = 1


    if len(plateDict['READS']) > 0:
        # Fix Barcode  ---------------------------------------------------------
        if 'BARCODE' in plateDict:
            if plateDict['BARCODE'] == 'Unable to read':
                plateDict.pop('BARCODE')
        
        # Set PlateID ---------------------------------------------------------
        if 'BARCODE' in plateDict:
            plateDict['PLATE_ID'] =  plateDict['BARCODE']
        else:   
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+xSheet
            else:
                plateDict['PLATE_ID'] = xSheet
        if len(plateDict['PLATE_ID']) > max_PlateID_len:
            logger.warning(f"WARNING : PlateID - length: {len(plateDict['PLATE_ID'])} / {max_PlateID_len}")
 
        # Setup TestPlate ---------------------------------------------------------
        djPlate = TestPlate.get(plateDict['PLATE_ID'],WellData=True,verbose=1)
        if djPlate is None:
            logger.info(f" New Plate {plateDict['PLATE_ID']}")
            djPlate = TestPlate.new(plateDict['PLATE_ID'],nWells,WellData=True)


        djPlate.reader = 'Tecan M1000 (iControl)'
        djPlate.experiment = ""
        djPlate.protocol = ""
        djPlate.test_operator = plateDict['USER'].split('\\')[-1]
        djPlate.test_date = plateDict['TEST_DATE']
        djPlate.n_reads = nWells

        # Fix ReadOut_Types ---------------------------------------------------------
        djPlate.n_readouts = len(plateDict['READS'])

        if "Absorbance" in plateDict['MODE'] :
            if djPlate.n_readouts > 1 :
                _mWaveLengths = plateDict['MEASUREMENT WAVELENGTH'].split(";")

                # Absorbance Resazurin OD (570-600)
                if '570' in plateDict['MEASUREMENT WAVELENGTH'] and '600' in plateDict['MEASUREMENT WAVELENGTH']:
                    _aR = _mWaveLengths.index('570')
                    _bR = _mWaveLengths.index('600')

                # Absorbance Readout_A - Readout_B
                else:
                    _aR = 0
                    _bR = 1

                # Set _readout_types for Wells  
                _readout_types = [f"OD{_mWaveLengths[_aR]}-{_mWaveLengths[_bR]}"]
                for w in _mWaveLengths:
                    _readout_types.append(f"OD{w}")

            else:
                # Absorbance OD (Wavelength)
                _readout_types = [f"OD{plateDict['MEASUREMENT WAVELENGTH']}"]     
            
        elif "Fluorescence" in plateDict['MODE'] :
            
            if "Bottom" in plateDict['MODE'] :
                _r = "Fb"
            elif "Top" in plateDict['MODE'] : 
                _r = "Ft"
            else:
                _r = 'F'
            _readout_types = [f"{_r}{plateDict['EXCITATION WAVELENGTH']}/{plateDict['EMISSION WAVELENGTH']}"] 

        djPlate.readout_type = _readout_types[0]

        # Well ReadOuts ---------------------------------------------------------
        for r in range(_nrow):
            for c in range(_ncol):
                _wellid = djPlate.well_id((r+1,c+1))
                if djPlate.n_readouts > 1:
                    _readouts = [plateDict['READS'][_aR][r][c] - plateDict['READS'][_bR][r][c]]
                    _readouts.append(plateDict['READS'][0][r][c])
                    _readouts.append(plateDict['READS'][1][r][c])
                else:
                    _readouts = [plateDict['READS'][0][r][c]]

                # Convert to DecimalField with 5 decimal points
                for i in range(len(_readouts)):
                    _readouts[i] = Decimal(_readouts[i]).quantize(Decimal("1.00000"))

                djPlate.wells[_wellid].readouts = _readouts
                djPlate.wells[_wellid].readout_types = _readout_types

        return(djPlate)
    else:
        return(None)
    
#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - BioTek HTC XLSX.SHEET
#--------------------------------------------------------------------------------------------------------------
def readPlate_Gen5_sheet(xSheet,xDF,prefix=None,):
    max_PlateID_len = 25

    plateDict = {}
    plateDict['READOUTS'] = []

    paramLst = ['Experiment File Path','Protocol File Path','Plate Number','Plate ID','Barcode',
                'Reader Type','Plate Type',
                'Measurement 1','Measurement 2',
#                'Excitation Wavelength','Emission Wavelength','Excitation Bandwidth','Emission Bandwidth',
#                'Gain','Number of Flashes','Flash Frequency','Integration Time','Lag Time',
#                'Settle Time'
                ]
    nWells = 0
    _flMatrix = 0
    _matrix = []
    _ncol = 0

    _nSheet = len(xDF)
    for r in range(_nSheet):
        #print(xDF[r])
        for k in paramLst:
            if k in str(xDF[0][r]):
                if str(xDF[1][r]) != "nan":
                    if k.upper() in plateDict:
                        plateDict[k.upper()] += ";"+str(xDF[1][r])
                    else:
                        plateDict[k.upper()] = str(xDF[1][r])

        if 'Reading Date' in str(xDF[0][r]):
            try:
                plateDict['TEST_DATE'] = datetime.strptime(str(xDF[1][r]), '%Y-%m-%d %H:%M:%S')
            except:
                plateDict['TEST_DATE'] = datetime.strptime(str(xDF[1][r]), '%Y-%d-%m %H:%M:%S')

        # MatrixRead - Row 
        if _flMatrix > 0 :
            _row = []
            if str(xDF[1][r]) in Plate.rowLabels:
                for c in range(_ncol):
                    _row.append(xDF[c+2][r])
                _matrix.append(_row)
            else:
                plateDict['READOUTS'].append(_matrix)
                _flMatrix = 0
                _matrix = []

        # MatrixRead - Init
        _rH = r+1
        if _rH > _nSheet-1:     # In case r+1 does not exists
            _rH = _nSheet-1

        if xDF[2][r] ==  1.0 and xDF[1][r+1] ==  'A':
            if xDF[13][r] == 12:
                nWells = 96
                _ncol = 12
                _nrow = 8
                _flMatrix = 1
            if xDF[25][r] == 24:
                nWells = 384
                _ncol = 24
                _nrow = 16
                _flMatrix = 1

    # In case matrix goes to last line
    if _flMatrix > 0 :
        plateDict['READOUTS'].append(_matrix)
        _flMatrix = 0
        _matrix = []

    if len(plateDict['READOUTS']) > 0:

        plateDict['MODE'] = "Absorbance"

        # Fix Barcode 
        if 'BARCODE' in plateDict:
            if plateDict['BARCODE'] == 'Unable to read':
                plateDict.pop('BARCODE')
        
        # Set PlateID
        if 'BARCODE' in plateDict:
            plateDict['PLATE_ID'] =  plateDict['BARCODE']
        else:   
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+xSheet
            else:
                plateDict['PLATE_ID'] = xSheet
        if len(plateDict['PLATE_ID']) > max_PlateID_len:
            print(f"WARNING : PlateID - length: {len(plateDict['PLATE_ID'])} / {max_PlateID_len}")

 
        # Setup TestPlate
        xPl = Plate(plateDict['PLATE_ID'],nWells,plateType='TestPlate')                                

        xPl.PlateData['PLATE_SIZE'] = f'{nWells}w'
        xPl.PlateData['TEST_DATE'] = plateDict['TEST_DATE']
        #xPl.PlateData['TEST_OPERATOR'] = plateDict['USER'].split('\\')[-1]
        xPl.PlateData['READER'] = 'BioTek HTX (Gen5)'
        xPl.PlateData['EXPERIMENT'] = os.path.split(plateDict['EXPERIMENT FILE PATH'])[1]
        xPl.PlateData['PROTOCOL']   = os.path.split(plateDict['PROTOCOL FILE PATH'])[1]
        xPl.PlateData['NREADS'] = len(plateDict['READOUTS'])
        #xPl.PlateData['PROTOCOL'] = 'OD600-Corning'
        xPl.PlateData['HAS_READOUT'] = 1

        # print(xPl.PlateData)
        # print(plateDict)
        # ReadOuts
        if xPl.PlateData['NREADS'] > 1:

            # Define calculation of multiple Readouts
            _aR = 0
            _bR = 1
            _calcType = 'Substraction'
            _mWaveLength = [plateDict['MEASUREMENT 1'],plateDict['MEASUREMENT 2']]

            if "Absorbance" in plateDict['MODE']:
                # fix OD570-600
                if "-".join(_mWaveLength) == '600-570':
                    _aR = 1
                    _bR = 0

            plateDict['MEASUREMENT'] = f"{_mWaveLength[_aR]}-{_mWaveLength[_bR]}"
            for r in range(_nrow):
                for c in range(_ncol):
                    xPl.set_WellProperty((r+1,c+1),'WELL_ID', xPl.map_WellID((r+1,c+1)))
                    xPl.set_WellProperty((r+1,c+1),'READOUT',plateDict['READOUTS'][_aR][r][c]-plateDict['READOUTS'][_bR][r][c])
                    xPl.set_WellProperty((r+1,c+1),'READOUTA',plateDict['READOUTS'][_aR][r][c])
                    xPl.set_WellProperty((r+1,c+1),'READOUTB',plateDict['READOUTS'][_bR][r][c])

        else:
            for r in range(_nrow):
                for c in range(_ncol):
                    xPl.set_WellProperty((r+1,c+1),'WELL_ID', xPl.map_WellID((r+1,c+1)))
                    xPl.set_WellProperty((r+1,c+1),'READOUT',plateDict['READOUTS'][0][r][c])
            plateDict['MEASUREMENT'] = plateDict['MEASUREMENT 1']

        # Fix ReadOut_ID 
        if "Absorbance" in plateDict['MODE'] :
            plateDict['READOUT_ID'] = f"OD{plateDict['MEASUREMENT']}"
        # elif plateDict['MODE'] == "Fluorescence Bottom Reading":
        #     plateDict['READOUT_ID'] = f"Fb{plateDict['EXCITATION WAVELENGTH']}/{plateDict['EMISSION WAVELENGTH']}" 
        # elif plateDict['MODE'] == "Fluorescence Top Reading":
        #     plateDict['READOUT_ID'] = f"Ft{plateDict['EXCITATION WAVELENGTH']}/{plateDict['EMISSION WAVELENGTH']}"
        xPl.PlateData['READOUT_ID'] = plateDict['READOUT_ID']

        #print(plateDict['READOUTS'])
        return(xPl)
    else:
        return(None)

#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - BMG CLARIOstar XLSX.SHEET
#--------------------------------------------------------------------------------------------------------------
def readPlate_BMG_sheet(xSheet,xDF,prefix=None,):
    max_PlateID_len = 25

    plateDict = {}
    plateDict['READOUTS'] = []

    paramLst = ['Experiment File Path','Protocol File Path','Plate Number','Plate ID','Barcode',
                'Reader Type','Plate Type',
                'Measurement 1','Measurement 2',
#                'Excitation Wavelength','Emission Wavelength','Excitation Bandwidth','Emission Bandwidth',
#                'Gain','Number of Flashes','Flash Frequency','Integration Time','Lag Time',
#                'Settle Time'
                ]
    nWells = 0
    _flMatrix = 0
    _matrix = []
    _ncol = 0

    _sDate = False
    _sTime = False
    
    _nSheet = len(xDF)
    for r in range(_nSheet):
        if 'Date' in str(xDF[0][r]):
            _sDate = xDF[0][r].replace('Date: ','')
        if 'Time' in str(xDF[0][r]):
            _sTime = xDF[0][r].replace('Time: ','')
        if _sDate and _sTime:
            plateDict['TEST_DATE'] = datetime.strptime(f"{_sDate} {_sTime}", '%d/%m/%Y %H:%M:%S %p')    

        if 'Test Name' in str(xDF[0][r]):
            plateDict['EXPERIMENT'] = xDF[0][r].replace('Test Name: ','')

        if 'ID2' in str(xDF[0][r]):
            if 'OD405' in xDF[0][r]:
                plateDict['PROTOCOL'] = 'OD405_Haemolysis.1'
                plateDict['READOUT_ID'] = 'OD405'

        # MatrixRead - Row 
        if _flMatrix > 0 :
            _row = []
            if str(xDF[0][r]) in Plate.rowLabels:
                for c in range(_ncol):
                    _row.append(xDF[c+1][r])
                _matrix.append(_row)
            else:
                plateDict['READOUTS'].append(_matrix)
                _flMatrix = 0
                _matrix = []

        # MatrixRead - Init
        _rH = r+1
        if _rH > _nSheet-1:     # In case r+1 does not exists
            _rH = _nSheet-1

        if xDF[1][r] ==  1.0 and xDF[0][r+1] ==  'A':
            if xDF[12][r] == 12:
                nWells = 96
                _ncol = 12
                _nrow = 8
                _flMatrix = 1
            if xDF[24][r] == 24:
                nWells = 384
                _ncol = 24
                _nrow = 16
                _flMatrix = 1

    # In case matrix goes to last line
    if _flMatrix > 0 :
        plateDict['READOUTS'].append(_matrix)
        _flMatrix = 0
        _matrix = []


    if len(plateDict['READOUTS']) > 0:
        # Fix Barcode 
        if 'BARCODE' in plateDict:
            if plateDict['BARCODE'] == 'Unable to read':
                plateDict.pop('BARCODE')
        
        # Set PlateID
        if 'BARCODE' in plateDict:
            plateDict['PLATE_ID'] =  plateDict['BARCODE']
        else:   
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+xSheet
            else:
                plateDict['PLATE_ID'] = xSheet
        if len(plateDict['PLATE_ID']) > max_PlateID_len:
            print(f"WARNING : PlateID - length: {len(plateDict['PLATE_ID'])} / {max_PlateID_len}")

        # Setup TestPlate
        xPl = Plate(plateDict['PLATE_ID'],nWells,plateType='TestPlate')                                

        xPl.PlateData['PLATE_SIZE'] = f'{nWells}w'
        xPl.PlateData['TEST_DATE'] = plateDict['TEST_DATE']
        #xPl.PlateData['TEST_OPERATOR'] = plateDict['USER'].split('\\')[-1]
        xPl.PlateData['READER'] = 'BMG CLARIOstar'
        xPl.PlateData['NREADS'] = len(plateDict['READOUTS'])
        xPl.PlateData['PROTOCOL'] = plateDict['PROTOCOL']
        xPl.PlateData['HAS_READOUT'] = 1

        # ReadOuts
        for r in range(_nrow):
            for c in range(_ncol):
                xPl.set_WellProperty((r+1,c+1),'WELL_ID', xPl.map_WellID((r+1,c+1)))
                xPl.set_WellProperty((r+1,c+1),'READOUT',plateDict['READOUTS'][0][r][c])
        xPl.PlateData['READOUT_ID'] = plateDict['READOUT_ID']

#        print(xPl.PlateData)


        #print(plateDict['READOUTS'])
        return(xPl)
    else:
        return(None)