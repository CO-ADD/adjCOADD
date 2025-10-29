
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
def multimodereader_xls(xlFile, prefix=None, as_is=False, **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    fXlsx = open(xlFile, "rb")
    xlWB = pd.ExcelFile(fXlsx)
    #xlWB = pd.ExcelFile(xlFile)
    
    dictPlates = {}
    for xSheet in xlWB.sheet_names:
        
        xDF = xlWB.parse(xSheet, header = None)

        if len(xDF)>0:
            sng_plate = {}
            djTP = None
            if xDF[0][0] == "Application: Tecan i-control":
                djTP,_status = read_iControl_xlsheet(xSheet,xDF,prefix=prefix)
                
            elif xDF[0][0] == "Experiment" or xDF[0][1] == "Experiment":
                djTP,_status = read_Gen5_sheet(xSheet,xDF,prefix=prefix)
            
            # elif "CLARIOstar" in xDF[0][3]:
            #     xPl = readPlate_BMG_sheet(xSheet,xDF,prefix=prefix)

            if djTP:
                djTP.input_file = os.path.split(xlFile)[1]

                sng_plate['plate_id'] = djTP.plate_id
                sng_plate['plate'] = djTP
                sng_plate['new'] = _status == 'New'
                
                if djTP.plate_id in dictPlates:
                    sng_plate['new'] = 'Duplicate'
                    _status = "Duplicate"
                else:
                    dictPlates[djTP.plate_id] = sng_plate
                    
                if 'SHEET' in djTP.plate_id.upper():
                    _status = "No Barcode"

                # Output Verbose/valLog
                if verbose>0:
                    logger.info(f"[{djTP.plate_id:25s}] - {djTP.reader}  {djTP.n_wells}w {djTP.readout_type} [{_status}]")
                if valLog:
                    if _status == 'New':
                        valLog.add_info("New Testplate",
                                        djTP.plate_id, 
                                        f"{djTP.reader}  {djTP.n_wells}w {djTP.readout_type}",
                                        "Select Upload")
                    elif _status == 'Exists':
                        valLog.add_warning("TestPlate Exists",
                                           djTP.plate_id,
                                           f"{djTP.reader}  {djTP.n_wells}w {djTP.readout_type}",
                                           "Select Overwrite")
                    elif _status == 'Duplicate':
                        valLog.add_error("Duplicate Testplate",
                                         djTP.plate_id,
                                         f"{djTP.reader}  {djTP.n_wells}w {djTP.readout_type}",
                                         "Correct PlateID in  Xlsx file")
                    elif _status == 'No Barcode':
                        valLog.add_error("No Barcode",
                                         djTP.plate_id,
                                         f"{djTP.reader}  {djTP.n_wells}w {djTP.readout_type}",
                                         "Correct PlateID in Xlsx file")
            else:
                if verbose>0:
                    logger.info(f"[{xSheet}] - Unknown PlateReader Format")
                if valLog:
                    valLog.add_warning("Unknown PlateReader Format",
                               f" Xls.Sheet: {xSheet}","",
                               "Check Xls.Sheet")                 
    fXlsx.close()
    return(list(dictPlates.values()))


#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - - Tecan iControl XLSX.SHEET
#--------------------------------------------------------------------------------------------------------------
def read_iControl_xlsheet(xSheet,xDF,prefix=None):
    max_PlateID_len = 25

    # Read XLS Sheet into plateDict ------------------------------------------------
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
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+plateDict['BARCODE']
            else:
                plateDict['PLATE_ID'] =  plateDict['BARCODE']
        else:   
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+xSheet
            else:
                plateDict['PLATE_ID'] = xSheet
        if len(plateDict['PLATE_ID']) > max_PlateID_len:
            logger.warning(f"WARNING : PlateID - length: {len(plateDict['PLATE_ID'])} / {max_PlateID_len}")
 

        # Set ReadOut_Types ---------------------------------------------------------    
        if "Absorbance" in plateDict['MODE'] :
            if len(plateDict['READS']) > 1 :
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

        # Setup TestPlate ---------------------------------------------------------
        _status = "Exists"
        djPlate = TestPlate.get(plateDict['PLATE_ID'],WellData=True,verbose=0)
        if djPlate is None:
            #logger.info(f" New Plate {plateDict['PLATE_ID']}")
            djPlate = TestPlate.new(plateDict['PLATE_ID'],nWells,WellData=True)
            _status = "New"

        djPlate.reader = 'Tecan M1000 (iControl)'
        djPlate.experiment = ""
        djPlate.protocol = ""
        djPlate.test_operator = plateDict['USER'].split('\\')[-1]
        djPlate.test_date = plateDict['TEST_DATE']
        djPlate.n_reads = nWells

        djPlate.n_readouts = len(_readout_types)
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
                #for i in range(len(_readouts)):
                #    _readouts[i] = Decimal(_readouts[i]).quantize(Decimal("1.00000"))

                djPlate.wells[_wellid].readouts = _readouts
                djPlate.wells[_wellid].readout_types = _readout_types

        return(djPlate,_status)
    else:
        return(None,False)

#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - Tecan iControl XML
#--------------------------------------------------------------------------------------------------------------
def read_iControl_xml(xmlFile,fileBarcode=None): 

    if os.path.exists(xmlFile):
        tree = ElementTree.parse(xmlFile)
        root = tree.getroot()
        Barcode = None
        for plate in root.iter('Plate'):
            for p in plate.iter():
                if p.tag == 'BC':
                    Barcode = p.text.strip()
                    if Barcode == '':
                        Barcode = None
                    #print(f"[{Barcode}]")
        for script in root.iter('Script'):
            for s in script.iter():
                if s.tag == '{tecan.at.schema.documents}ReadingLabel':
                    ReadLabel = s.attrib['name']
                if Barcode is None:
                    if s.tag == '{tecan.at.schema.documents}Barcode':
                        if 'name' in s.attrib: 
                            Barcode = s.attrib['name']
                #print(s.tag)
        for section in root.iter('Section'):
            #print(section.attrib['Time_End'])
            TestDate = pd.to_datetime(section.attrib['Time_End'])
            #print(datetime.fromisoformat(section.attrib['Time_End']))
            for parameter in section.iter('Parameter'):
                #print(parameter.attrib)
                if parameter.attrib['Name'] == 'Mode':
                    Mode = parameter.attrib['Value']
                    #print()
                if parameter.attrib['Name'] == 'Wavelength':
                    Wavelength = parameter.attrib['Value']
                    WavelengthUnit = parameter.attrib['Unit']
            rWell = {}
            for well in section.iter('Well'):
                #print(well.attrib['Pos'])
                #print(well[0].text)
                rWell[well.attrib['Pos']] = well[0].text
            nWells = len(rWell)
            PlateSize = f"{nWells}w"
            if Barcode is None:
                #print(f"No Barcode -> {fileBarcode}")
                if fileBarcode:
                    Barcode = fileBarcode
                else:
                    Barcode,_ = os.path.splitext(os.path.split(xmlFile)[1])

        xPl = Plate(Barcode,nWells)

        readMode = ""
        if Mode == "Absorbance":
            readMode = "OD"
        if Wavelength:
            readMode += Wavelength
        xPl.PlateData['PLATE_SIZE'] = PlateSize
        xPl.PlateData['TEST_DATE'] = TestDate
        xPl.PlateData['READER'] = 'Tecan M1000 (iControl)'
        xPl.PlateData['NREADS'] = 1
        xPl.PlateData['PROTOCOL'] = 'OD600-Corning'
        xPl.PlateData['READOUT_ID'] = readMode
        xPl.PlateData['INPUTFILE'] = os.path.split(xmlFile)[1]
        xPl.PlateData['HAS_READOUT'] = 1

        for w in rWell:
            xPl.set_WellProperty(w,'READOUT',rWell[w])
        print(f" [XML-iControl] {Barcode} {PlateSize} [{xmlFile}]")
        return(xPl)
    return(None)
    
#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - BioTek HTC XLSX.SHEET
#--------------------------------------------------------------------------------------------------------------
def read_Gen5_sheet(xSheet,xDF,prefix=None,):
    max_PlateID_len = 25

    # Read XLS Sheet into plateDict ------------------------------------------------
    plateDict = {}
    plateDict['READS'] = []

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
            if str(xDF[1][r]) in TestPlate.ROW_LABELS:
                for c in range(_ncol):
                    _row.append(xDF[c+2][r])
                _matrix.append(_row)
            else:
                plateDict['READS'].append(_matrix)
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
        plateDict['READS'].append(_matrix)
        _flMatrix = 0
        _matrix = []

    if len(plateDict['READS']) > 0:
        # Fix Barcode  ---------------------------------------------------------
        if 'BARCODE' in plateDict:
            if plateDict['BARCODE'] == 'Unable to read':
                plateDict.pop('BARCODE')
        
        # Set PlateID ---------------------------------------------------------
        #print(plateDict)
        if 'BARCODE' in plateDict:
            plateDict['PLATE_ID'] =  plateDict['BARCODE']
        elif 'PLATE ID' in plateDict:
            plateDict['PLATE_ID'] =  plateDict['PLATE ID']
        else:   
            if prefix:
                plateDict['PLATE_ID'] = prefix+"_"+xSheet
            else:
                plateDict['PLATE_ID'] = xSheet
        if len(plateDict['PLATE_ID']) > max_PlateID_len:
            logger.warning(f"WARNING : PlateID - length: {len(plateDict['PLATE_ID'])} / {max_PlateID_len}")
 
        # Set ReadOut_Types ---------------------------------------------------------    
        plateDict['MODE'] = "Absorbance"

        if "Absorbance" in plateDict['MODE'] :
            if len(plateDict['READS']) > 1 :
                _mWaveLengths = plateDict['MEASUREMENT 1'],plateDict['MEASUREMENT 2']

                # Absorbance Resazurin OD (570-600)
                if '570' in _mWaveLengths and '600' in _mWaveLengths:
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
                _readout_types = [f"OD{plateDict['MEASUREMENT 1']}"]     


        # Setup TestPlate ---------------------------------------------------------
        _status = "Exists"
        djPlate = TestPlate.get(plateDict['PLATE_ID'],WellData=True,verbose=0)
        if djPlate is None:
            #logger.info(f" New Plate {plateDict['PLATE_ID']}")
            djPlate = TestPlate.new(plateDict['PLATE_ID'],nWells,WellData=True)
            _status = "New"

        djPlate.reader = 'BioTek HTX (Gen5)'
        djPlate.experiment = os.path.split(plateDict['EXPERIMENT FILE PATH'])[1]
        djPlate.protocol = os.path.split(plateDict['PROTOCOL FILE PATH'])[1]
        #djPlate.test_operator = plateDict['USER'].split('\\')[-1]
        djPlate.test_date = plateDict['TEST_DATE']
        djPlate.n_reads = nWells

        djPlate.n_readouts = len(_readout_types)
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
                #for i in range(len(_readouts)):
                #    _readouts[i] = Decimal(_readouts[i]).quantize(Decimal("1.00000"))

                djPlate.wells[_wellid].readouts = _readouts
                djPlate.wells[_wellid].readout_types = _readout_types

        return(djPlate,_status)
    else:
        return(None,False)
        
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

#--------------------------------------------------------------------------------------------------------------
# Plate - Read TestPlate Readout - BioTek Epoch 
#--------------------------------------------------------------------------------------------------------------
def readPlate_Epoch_txt(pltFile,fileBarcode=None): 

    if os.path.exists(pltFile):
        with open(pltFile,"r") as file:
            txtData = list(csv.reader(file, delimiter= '\t'))

        plSize = 0
        plDate = None
        plProtocol = None

        # Get Plate Size
        if len(txtData[1]) > 24:
            if txtData[1][24] == '24' and txtData[17][0] == 'P':
                plSize = 384
        else:
            if txtData[1][12] == '12' and txtData[9][0] == 'H':
                plSize = 96

        # Get last Datetime
        flDate = False
        for r in txtData:
            if len(r)>0:
                if flDate:
                    plDate = pd.to_datetime(r[0])
                    #plDate = datetime.strptime(r[0],"%d/%m/%Y %H:%M:%S %p")
                if r[0] == 'Date':
                    flDate=True

        # Get last Protocol
        for r in txtData:
            if len(r)>2:
                if 'Protocol' in r[2]:
                    res = re.findall(r'([A-Za-z0-9\s]+).prt',r[2])
                    plProtocol = f"{res[0]}.prt" 

        readMode = "OD" + txtData[0][0]
        if fileBarcode:
            Barcode = fileBarcode
        else:
            Barcode = os.path.splitext(os.path.split(pltFile)[1])[0]

        if plSize > 0:
            xPl = Plate(Barcode,plSize)
            xPl.PlateData['PLATE_SIZE'] = f"{plSize}w"
            xPl.PlateData['TEST_DATE'] = plDate
            xPl.PlateData['READER'] = 'Biotek Epoch (Gen5)'
            xPl.PlateData['NREADS'] = 1
            xPl.PlateData['PROTOCOL'] = plProtocol
            xPl.PlateData['READOUT_ID'] = readMode
            xPl.PlateData['INPUTFILE'] = os.path.split(pltFile)[1]
            xPl.PlateData['HAS_READOUT'] = 1

            readRow = 1
            readCol = 0
            for r in range(1,xPl.rows+1):
                for c in range(1,xPl.columns+1):
                    xPl.set_WellProperty((r,c),'READOUT',txtData[readRow+r][readCol+c])

            print(f" [TXT-Epoch] {Barcode} {plSize}w [{pltFile}]")
            return(xPl)
    return(None)