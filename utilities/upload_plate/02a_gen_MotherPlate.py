import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from functools import reduce
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "genMotherPlate"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
#    format="[%(name)-20s] %(message)s ",
    format="%(message)s",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------


def gen_MotherPlate_HCR(MP_Dict,Racks,plate_size=384):
    Q = {'A1':(0,0),'B1':(1,0),'A2':(0,1),'B2':(1,1)}
    Properties=['BARCODE','COMPOUND_ID','PROJECT_ID','COMPOUND_CODE',
                'CONC','CONC_UNIT','SOLVENT','VOLUME','VOLUME_UNIT',
                'FULL_MW','FULL_MF']
    
    #print(MP_Dict)
    #print(Racks)
    _plate_id = MP_Dict['MOTHERPLATEID']
    #_mp = Plate(_plate_id,plate_size,plateData=None,plateType='MasterPlate')

    _rack_start = 0
    if MP_Dict['RACKROW'] == 'I':
        _rack_start = 4
    
    for _rack_row in range(1,5):
        for _rcol in range(1,13):
            _rrow = _rack_row + _rack_start
            _mrow = (int((_rack_row - 1)/2) * 8) + 1
            _mcol = ((_rcol - 1)*2) + ((_rrow - 1) % 2) + 1
            
            _rpos = Racks[MP_Dict['RACKID']].map_pos1D((_rrow,_rcol))
            _rwell = Racks[MP_Dict['RACKID']].WellData[_rpos]

            _mpos = _mp.map_pos1D((_mrow,_mcol))

            if _rwell['BARCODE'] != 'NO READ':
                #print(f"[{_rack_row}] {_rrow} {_rcol} -> {_mrow} {_mcol}")
                #print(f"[{_rack_row}] {Racks[MP_Dict['RACKID']].map_WellID((_rrow,_rcol))} -> {_mp.map_WellID((_mrow,_mcol))}")
                
                for col in Properties:
                    _mp.set_WellProperty(_mpos,col,_rwell[col])
                _mp.set_WellProperty(_mpos,MP_Dict['RACKID'],'RackID')    
                _mp.set_WellProperty(_mpos,Racks[MP_Dict['RACKID']].map_WellID((_rrow,_rcol)),'WellID')    
    return(_mp)

def gen_MotherPlate_PSR(MP_Dict,Racks,plate_size=384):
    Q = {'A1':(0,0),'B1':(1,0),'A2':(0,1),'B2':(1,1)}
    Properties=['BARCODE','COMPOUND_ID','PROJECT_ID','COMPOUND_CODE',
                'CONC','CONC_UNIT','SOLVENT','VOLUME','VOLUME_UNIT',
                'FULL_MW','FULL_MF']
    
    _plate_id = MP_Dict['MOTHERPLATEID']
    _mp = Plate(_plate_id,plate_size,plateData=None,plateType='MasterPlate')

    _setid = 1
    for qq in Q.keys():
        if MP_Dict[qq]:
            _rack = Racks[MP_Dict[qq]]
            #map96to384(_rack,_mp,qq,Properties=['BARCODE','COMPOUND_ID','PROJECT_ID','COMPOUND_CODE'])
            
            for pos,well in _rack.WellData.items():
                if well['BARCODE'] != 'NO READ':

                    # Map 96 to 384 by Quadrants
                    _r,_c = _rack.map_pos2D(well['WELL_ID']) 
                    _qr = ((_r - 1) * 2) + 1 + Q[qq][0]
                    _qc = ((_c - 1) * 2) + 1 + Q[qq][1]
                    _mpos = _mp.map_WellID((_qr,_qc))

                    # Set the MotherWell Properties
                    for col in Properties:
                        _mp.set_WellProperty(_mpos,col,well[col])
                    _mp.set_WellProperty(_mpos,'SET_ID',str(_setid))
                    _mp.set_WellProperty(_mpos,MP_Dict['RACKID'],'RackID')    
                    _mp.set_WellProperty(_mpos,Racks[MP_Dict['RACKID']].map_WellID((_rrow,_rcol)),'WellID')    
    _setid += 1
    return(_mp)


# --------------------------------------------------------------------------------
def read_hcprep_prepsheet_xls(xlFile, SheetName='HCPrep', prefix=None, as_is=False, **kwargs):
# --------------------------------------------------------------------------------
    xlWB = pd.ExcelFile(xlFile)
    xDF = xlWB.parse(SheetName)
    xDF.columns = [c.lower() for c in xDF.columns]

    print(xDF)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from applib.plate.masterplates import read_motherplate_prepsheet_xls
    from dscreen.models import Screen_Run
    from dplate.models import MasterPlate, MasterWell
    from dsample.models import COADD_Compound, ABase_Compound_Batch
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # MotherPlate XLSX -------------------------------------------------------------
    verbose = 0

    if prgArgs.rackdir and prgArgs.runid and prgArgs.excelfile:
        
        djRun = Screen_Run.get(prgArgs.runid)
        if djRun is None:
            djRun = Screen_Run()
            djRun.run_id = prgArgs.runid
            new_runid = True
            print(f" [{prgArgs.table}] RunID {prgArgs.runid} Does NOT Exist ")

        RackDir = prgArgs.rackdir
        rack_files = [f for f in os.listdir(RackDir) if f.endswith(".csv")]

        # Dataframe of all the Rack Files
        Racks = {}
        for rack_file in rack_files:
            _n_tubes = 0
            _n_cmpds = 0
            _rack_df = pd.read_csv(os.path.join(RackDir,rack_file))
            _rack_df.columns =  [c.upper() for c in _rack_df.columns]
            _rack_id = str(_rack_df['RACKID'].unique()[0])

            #_rack_df=_rack_df.apply(apply_get_barcode, axis=1)
            #print(f" [Rack] {_rack_id} from {rack_file}")

            _rack = MasterPlate.new(_rack_id,96,'Storage',WellData=True)
            for idx,row in _rack_df.iterrows():
                if row['BARCODE'] != 'NO READ':
                    _n_tubes += 1
                    _well = _rack.get_well(row['WELLID'])
                    _well.barcode = row['BARCODE']

                    _tube = MasterWell.get(None,None,row['BARCODE'])
                    if _tube:
                        _n_cmpds += 1
                        _well.cmpbatch_id = _tube.cmpbatch_id
                    else:
                        _well.cmpbatch_id = None

            Racks[_rack_id] = _rack

            print(f" [{prgArgs.table}] Reading {rack_file} [Barcodes: {_n_tubes}, Cmpds: {_n_cmpds}] ")
        
        # for _rack_id in Racks:
        #     print(Racks[_rack_id])

        # == HC-Prep =======================================================================================
        if prgArgs.table == 'HCPrep':
            MP_Prep_Xlsx = 'MotherPlate_from_HCPrep.xlsx'
    
            xlWB = pd.ExcelFile(prgArgs.excelfile)
            HCPrep = xlWB.parse('HCPrep')
            HCPrep.columns = [c.upper() for c in HCPrep.columns]

            # For each MotherPlate <- Rack from row A or E 
            MPs = {}
            for idx,row in HCPrep.iterrows():
                if row['MOTHERPLATEID'] not in MPs:
                    MPs[row['MOTHERPLATEID']]= MasterPlate.new(row['MOTHERPLATEID'],384,'Mother',WellData=True)

                _rack_start = -1
                if row['RACKROW'] == 'A':
                    _rack_start = 0
                elif row['RACKROW'] == 'E':
                    _rack_start = 4

                if _rack_start >= 0:
                    for _rack_row in range(1,5):
                        for _rcol in range(1,13):
                            # 'A' : 1st Row-> Left 1st Row, 2nd-> Right 1st Row
                            # 'I' : 3rd Row-> Left 9st Row, 4th-> Right 9st Row
                            _rrow = _rack_row + _rack_start
                            _mrow = (int((_rack_row - 1)/2) * 8) + 1
                            _mcol = ((_rcol - 1)*2) + ((_rrow - 1) % 2) + 1

                            _tube = Racks[row['RACKID']].get_well((_rrow,_rcol))
                            _mp_well = MPs[row['MOTHERPLATEID']].get_well((_mrow,_mcol))

                            _mp_well.barcode = _tube.barcode
                            _mp_well.cmpbatch_id = _tube.cmpbatch_id

                            if verbose>0:
                                print(f" {row['RACKID']} {_tube.well_id} -> {row['MOTHERPLATEID']} {_mp_well.well_id} [{_mp_well.barcode} {_mp_well.cmpbatch_id}] ")

        # == PS-Prep =======================================================================================

        if prgArgs.table == 'PSPrep':
            Q = {'A1':(0,0),'B1':(1,0),'A2':(0,1),'B2':(1,1)}

            MP_Prep_Xlsx = 'MotherPlate_from_PSPrep.xlsx'

            xlWB = pd.ExcelFile(prgArgs.excelfile)
            PSPrep = xlWB.parse('PSPrep')
            PSPrep.columns = [c.upper() for c in PSPrep.columns]

            # For each MotherPlate <- Racks A1, B1, A2, B2 
            MPs = {}
            for idx,row in PSPrep.iterrows():
                if row['MOTHERPLATEID'] not in MPs:
                    MPs[row['MOTHERPLATEID']]= MasterPlate.new(row['MOTHERPLATEID'],384,'Mother',WellData=True)

                _setid = 1
                for qq in Q.keys():
                    if row[qq]:
                        _rack = Racks[row[qq]]

                        for w in _rack.wells:
                            if _rack.wells[w].barcode:

                                # Map 96 to 384 by Quadrants
                                _r,_c = _rack.map_pos2D(w) 
                                _qr = ((_r - 1) * 2) + 1 + Q[qq][0]
                                _qc = ((_c - 1) * 2) + 1 + Q[qq][1]

                                _tube = _rack.get_well(w)
                                _mp_well = MPs[row['MOTHERPLATEID']].get_well((_qr,_qc))

                                _mp_well.barcode = _tube.barcode
                                _mp_well.cmpbatch_id = _tube.cmpbatch_id

                                # _mpos = _mp.map_WellID((_qr,_qc))

                                # # Set the MotherWell Properties
                                # for col in Properties:
                                #     _mp.set_WellProperty(_mpos,col,well[col])
                                # _mp.set_WellProperty(_mpos,'SET_ID',str(_setid))
                                # _mp.set_WellProperty(_mpos,MP_Dict['RACKID'],'RackID')    
                                # _mp.set_WellProperty(_mpos,Racks[MP_Dict['RACKID']].map_WellID((_rrow,_rcol)),'WellID')    
                _setid += 1

        # == Create MotherPlate output - to be copied into [MotherPlate] ===============================
        MP_Wells = []
        for _mp in MPs:
            for w in MPs[_mp].wells:
                if MPs[_mp].wells[w].barcode:
                    _mp_well_dict = {'MotherPlate_ID':_mp,
                                        'MotherWell_ID':w,
                                        'Plating':'IMB',
                                    }
                    
                    if MPs[_mp].wells[w].cmpbatch_id:
                        _mp_well_dict['CompoundID'] = str(MPs[_mp].wells[w].cmpbatch_id)
                        if MPs[_mp].wells[w].cmpbatch_id.batch_source == 'COADD':
                            _cmp = COADD_Compound.get(MPs[_mp].wells[w].cmpbatch_id)
                            _mp_well_dict['CompoundName'] = _cmp.compound_code
                            _mp_well_dict['ProjectID'] = str(_cmp.project_id)

                        elif MPs[_mp].wells[w].cmpbatch_id.batch_source == 'COADD':
                            _cmp = ABase_Compound_Batch.get(MPs[_mp].wells[w].cmpbatch_id)
                            _mp_well_dict['CompoundName'] = _cmp.compound_code
                            _mp_well_dict['ProjectID'] = str(_cmp.project_id)
                    else:
                        _mp_well_dict['CompoundID'] = "BARCODE NOT FOUND"
                        _mp_well_dict['CompoundName'] = ""
                        _mp_well_dict['ProjectID'] = ""

                    _mp_well_dict['Barcode']= MPs[_mp].wells[w].barcode
                    MP_Wells.append(_mp_well_dict)

        dfMP = pd.DataFrame(MP_Wells)
        print(f" [{prgArgs.table}] Write {MP_Prep_Xlsx} [MP Wells: {len(dfMP)}] ")
        dfMP.to_excel(MP_Prep_Xlsx)



#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t","--table",default=None,required=True, dest="table", action='store', help="Table to upload [HCPrep/PSPrep]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")

#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")

    prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
    prgParser.add_argument("-d","--rackdir",default=None,required=True, dest="rackdir", action='store', help="Directory for Rack files")
    prgParser.add_argument("-e","--excel",default=None,required=True, dest="excelfile", action='store', help="Excel File")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    from zDjango.djUtils import init_django_dir

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)

#==============================================================================
