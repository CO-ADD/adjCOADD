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

    if _prepSheets[PREP_SHEET] is not None:
        xDF = _prepSheets[PREP_SHEET]

        # Create/Get MasterPlates
        for _mpid in xDF['masterplate'].unique():
            _mp_new = True
            _rackid = str(_mpid)
            lstMP[_rackid] = {'plate_id':_rackid}

            _MP = MasterPlate.get(_rackid, WellData=True, verbose=0)
            if _MP is None:
                nWells=384
                _MP = MasterPlate.new(_rackid, nWells, PlateType='Storage', WellData=True)
                valLog.add_info('New Rack', _rackid,f'New Storage TubeRack',"Select Upload")
                _mp_new = True
            else:
                #_MP.load_wells(WellModel=MasterWell)
                _mp_new = False

            # for w in  _MP.wells:
            #     print(f"{_MP.wells[w]}")       

            lstMP[_rackid]['plate'] = _MP
            lstMP[_rackid]['new'] = _mp_new

        # For each Compound
        for idx,row in xDF.iterrows():
            
            _mw_status = "New"
            validStatus = True

            # Check if Barcode exists ----------------
            _barcode_status = "New"
            if row['barcode']:
                # In case of numeric only Barcode's
                row['barcode'] = str(row['barcode'])

                if MasterWell.exists(None,None,row['barcode']):
                    _barcode_status = "Exists"
                    valLog.add_error('Barcode exists', row['barcode'],f"Existing Barcode for {row['compound_id']}")

            #Check if Compound exists
            _cmpd_status = 'Exists'
            djCmp = Compound_Batch.get(row['compound_id'])
            if not djCmp:
                valLog.add_error('Wrong Compound_ID', row['compound_id'],
                    f"Compound (Batch) not found","Upload Compounds")
                _cmpd_status = "Missing"

            _mw_status = "New"
            _rackid = str(row['masterplate'])
            _wellid = row['masterwell']

            # Check if MasterWell has CmpdBatch_ID or Barcode
            _MP = lstMP[_rackid]['plate']
            _Well = _MP.get_well(_wellid)

            if _Well.barcode:
                valLog.add_error('Existing Rack/Pos', f"{_rackid}:{_wellid}",
                                        f"Rack has Tube in this position: {_Well.barcode} [{_Well.cmpbatch_id}]","Relocate Tube/Barcode first")
                _mw_status = 'Exists'
            elif _Well.cmpbatch_id:
                valLog.add_error('Existing Rack/Pos', f"{_rackid}:{_wellid}",
                                        f"Rack has Compound in this position: [{_Well.cmpbatch_id}]","Check MasterPlate/Well ID's")
                _mw_status = 'Exists'

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

                set_model_dicts(_Well,row,['conc_unit','amount_unit','solvent_conc_unit'])
                _Well.barcode = row['barcode']
                _Well.conc = Decimal(row['conc'])
                _Well.amount = Decimal(row['amount'])
                _Well.volume = Decimal(row['volume'])
                _Well.solvent = row['solvent']
                _Well.solvent_conc = Decimal(row['solvent_conc'])

                _Well.cmpbatch_lst = [row['compound_id']]
                _Well.cmpbatch_id = djCmp
                _Well.n_cmpbatches = 1
 
                # validDict = _Well.validate_model()       
                # if validDict:
                #     validStatus = False
                #     for c in validDict:
                #         print(f" [Upload_StockPrep] validDict Well {_rackid}:{_wellid} {c} ")
 
        #DEBUG output
        # for _mpid in lstMP:
        #     print(f" {_mpid} {repr(lstMP[_mpid]['plate'])}")

        #     for w in  lstMP[_mpid]['plate'].wells:
        #         print(f"{lstMP[_mpid]['plate'].wells[w]}")       

    else:
        if verbose>0:
            logger.error(f"{PREP_SHEET} not found in {xlFile}")
        if valLog:
            valLog.add_error("Wrong StockPrep file",
                            f"Sheet: {PREP_SHEET}", 
                            "SheetName not in StockPrep.XLS",
                            "Correct SheetName")        
    return(lstMP)
