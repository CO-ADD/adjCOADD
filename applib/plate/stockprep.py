import os
import pandas as pd
from decimal import Decimal

from dplate.models import MasterPlate, MasterWell
from dsample.models import Compound_Batch
from applib.data.set_fielddata import set_model_dicts

from django.conf import settings
import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
def get_StockPrep_xlsx(xlsFile, Sheets=[], FillNA='-', **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    CmpdPrep_Sheets = {
        'Samples': None,
        'StockPrep': None,
        }

    if os.path.isfile(xlsFile):
        fXlsx = open(xlsFile, "rb")
        xls = pd.ExcelFile(fXlsx)
        #xls = pd.ExcelFile(xlsFile)

        if Sheets is None:
            sheets = [CmpdPrep_Sheets.keys()]
        for key in Sheets:
            if key in CmpdPrep_Sheets:
                CmpdPrep_Sheets[key] = pd.read_excel(xls, key)
                CmpdPrep_Sheets[key].columns = [c.lower() for c in CmpdPrep_Sheets[key].columns]
                if FillNA:
                    CmpdPrep_Sheets[key] = CmpdPrep_Sheets[key].fillna(FillNA)
            else:
                if valLog:
                    valLog.add_error('Missing Sheet',key,f"XLSX {os.path.basename(xlsFile)}",f"Correct XLSX Sheets {list(CmpdPrep_Sheets)}")
        fXlsx.close()
        
    return(CmpdPrep_Sheets)

# --------------------------------------------------------------------------------
def read_Stock_Prepsheet_XLS(xlFile, prefix=None, **kwargs):
# --------------------------------------------------------------------------------

    PREP_SHEET = 'StockPrep'

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    lstMP = {}

    _prepSheets = get_StockPrep_xlsx(xlFile,Sheets=[PREP_SHEET],FillNA='-') 
    # if FileList:
    #     nFiles = len(FileList)
    # else:
    #     nFiles = 0
    # nUploads = 0

    # djPrj = Project.get(ProjectID)
    # valLog = Validation_Log("Upload_CmpdPrep")

             
    if _prepSheets[PREP_SHEET] is not None:
        xDF = _prepSheets[PREP_SHEET]

        for idx,row in xDF.iterrows():
            validStatus = True
            if row['barcode']:
                # In case of numeric only Barcode's
                row['barcode'] = str(row['barcode'])

                if MasterWell.exists(None,None,row['barcode']):
                    _barcode_status = "Exists"
                    valLog.add_error('Barcode exists', row['barcode'],f"Existing Barcode for {row['compound_id']}")
                else:
                    _barcode_status = "New"

            if row['masterplate'] and row['masterwell']:
                # In case of numeric only MasterPlate ID's
                row['masterplate'] = str(row['masterplate'])
                _mp_new = False
                _mw_status = "Exists"

                # Check MasterPlate
                if row['masterplate'] not in lstMP:
                    djMP = MasterPlate.get(row['masterplate'], WellData=True, verbose=0)
                    if djMP is None:
                        #logger.info(f" New Plate {mpid}")
                        nWells=384
                        djMP = MasterPlate.new(row['masterplate'], nWells, PlateType='Stock', WellData=True)
                        valLog.add_info('New Rack', row['masterplate'],f'New Storage TubeRack',"Select Upload")
                        _mp_new = True
                        _mw_status = "New"
                    else:
                        djMP.load_wells(WellModel=MasterWell)
                        djWell = djMP.get_well(row['masterwell'])

                    lstMP[row['masterplate']] = {'plate_id':row['masterplate'], 'plate': djMP, 'new':_mp_new}

                # Check MasterWell -- 
                djWell = lstMP[row['masterplate']]['plate'].get_well(row['masterwell'])
                if not lstMP[row['masterplate']]['new']:
                    if djWell.barcode:
                        valLog.add_error('Existing Rack/Pos', f"{row['masterplate']}:{row['masterwell']}",
                                            f"Rack has Tube in this position with {djWell.barcode}","Relocate Tube/Barcode first")
                    else:
                        valLog.add_warning('Existing Rack/Pos', f"{row['masterplate']}:{row['masterwell']}",
                                            f"Empty Position will be used for new Tube","Select Overwrite")
                        _mw_status = "Empty"

                #Check Compound
                djCmp = Compound_Batch.get(row['compound_id'])
                if not djCmp:
                    valLog.add_error('Wrong Compound_ID', row['compound_id'],
                        f"Compound (Batch) not found","Upload Compounds")
                    _mw_status = "No Compound"

                if _mw_status in ['New','Empty']:
                    if 'solvent_conc' not in row:
                        row['solvent_conc'] = 100
                    if 'solvent_conc_unit' not in row:
                        row['solvent_conc_unit'] = 'pct'

                    if 'amount_unit' not in row:
                        row['amount_unit'] = 'mg'
                    if 'volume_unit' not in row:
                        row['volume_unit'] = 'uL'
                    if 'conc_unit' not in row:
                        row['conc_unit'] = 'mg/mL'

                    set_model_dicts(djWell,row,['conc_unit','amount_unit','solvent_conc_unit'])
                    djWell.barcode = row['barcode']
                    djWell.conc = Decimal(row['conc'])
                    djWell.amount = Decimal(row['amount'])
                    djWell.volume = Decimal(row['volume'])
                    djWell.solvent = row['solvent']
                    djWell.solvent_conc = Decimal(row['solvent_conc'])

                    djWell.cmpbatch_lst = [row['compound_id']]
                    djWell.cmpbatch_id = djCmp
                    djWell.n_cmpbatches = 1

                    if settings.DEBUG:
                        print(f" {djWell.barcode} {djWell.cmpbatch_lst} {djWell}")
                    #     print(f" {djWell.conc} {djWell.amount} {djWell.volume} {djWell.solvent_conc}")

                    validDict = djWell.validate_model()
                    if validDict:
                        validStatus = False
                        for c in validDict:
                            print(f" [Upload_StockPre] validDict: {c} ")
    else:
        if verbose>0:
            logger.error(f"{PREP_SHEET} not found in {xlFile}")
        if valLog:
            valLog.add_error("Wrong StockPrep file",
                            f"Sheet: {PREP_SHEET}", 
                            "SheetName not in StockPrep.XLS",
                            "Correct SheetName")        
    return(lstMP)
