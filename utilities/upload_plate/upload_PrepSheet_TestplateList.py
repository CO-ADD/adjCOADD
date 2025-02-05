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
logName = "UploadTestPlateList"
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
def get_PlatePrep_xlsx(xlsFile, sheets=[]):

    PlatePrep_Sheets = {
        'TestPlateList': None,
        'MotherPlates': None,
        'Assays': None,
        'PSPrep': None,
        'HCPrep': None,
        }

    if os.path.isfile(xlsFile):
        xls = pd.ExcelFile(xlsFile)

        if sheets is None:
            sheets = [PlatePrep_Sheets.keys()]
        for key in sheets:
            if key in PlatePrep_Sheets:
                PlatePrep_Sheets[key] = pd.read_excel(xls, key)
                PlatePrep_Sheets[key].columns = [c.lower() for c in PlatePrep_Sheets[key].columns]
                
    return(PlatePrep_Sheets)


#-----------------------------------------------------------------------------
def load_PrepSheet_TestPlateList(db,xlPrepSheet,upload=False):
    dfTPl = oraCastDB.read_ExcelSheet(xlPrepSheet,"TestPlateList")
    dfTPl.columns = [c.upper() for c in dfTPl.columns]
    print(dfTPl)

    sCol = {"ASSAY_ID":"ASSAYTYPE_ID",
            "RESULT_TYPE":"RESULT_TYPE",
            "MOTHERPLATE_ID":"MOTHERPLATE_ID",
            "MOTHERPLATE2_ID":"MOTHERPLATE2_ID",
            "PLATING":"PLATING",
            "LABWARE":"LABWARE_ID",
            "TEST_STRAIN":"TEST_STRAIN",
            "TEST_MEDIA":"MEDIA_ID",
            "TEST_DYE":"TEST_DYE",
            "TEST_ADDITIVE":"TEST_ADDITIVE",
            "LAYOUT":"LAYOUT_CONTROL",
            "PROCESSING":"PROCESSING",
            "ISSUES":"ISSUES",
            "SYN_COMPOUNDS_AB":"SYN_COMPOUNDS_A",
            "SYN_COMPOUNDS_POT":"SYN_COMPOUNDS_B"
        }


    if len(dfTPl) > 0:
        for idx, row in dfTPl.iterrows():      
            fUpload = False
            nCnt = oraCastDB.check_PlateID_exists(db,"TestPlate",row['TESTPLATE_ID'])
            if nCnt == 1:
                fUpload = True
            else:
                fUpload = False
                logger.error(f"[CastDB] No TestPlate Found [{row['TESTPLATE_ID']}]")
            if pd.notnull(row['MOTHERPLATE_ID']):
                m1Cnt = oraCastDB.check_PlateID_exists(db,"MasterPlate",row['MOTHERPLATE_ID'])
                if m1Cnt == 1:
                    fUpload = True
                else:
                    fUpload = False
                    logger.error(f"[CastDB] No MotherPlate Found [{row['MOTHERPLATE_ID']}]")
            if pd.notnull(row['MOTHERPLATE2_ID']):
                if row['MOTHERPLATE2_ID']:
                    m2Cnt = oraCastDB.check_PlateID_exists(db,"MasterPlate",row['MOTHERPLATE2_ID'])
                    if m2Cnt == 1:
                        fUpload = True
                    else:
                        fUpload = False
                        logger.error(f"[CastDB] No MotherPlate2 Found [dlTPl[p]['MOTHERPLATE_ID']]")


            if fUpload and upload:
                cTable = "TestPlate"
                sWhere = f" Plate_ID = '{row['TESTPLATE_ID']}' "
                sDict = set_dictFields(row,sCol)
                sSql = db.gen_UpdateSQL(cTable,sDict,sWhere,bindvars=True)

                logger.info(f"[CastDB] Updating {row['TESTPLATE_ID']} PlateData ")
                db.exec(sSql,sDict,commit=True)
            else:
                logger.info(f"[CastDB] NO Updating for {row['TESTPLATE_ID']} ")

    iCol = "TestPlate_ID"
#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    django.setup()

    from dscreen.models import Assay
    from dplate.models import Labware, TestPlate, TestWell
    from dorganism.models import Organism, Organism_Batch
    from dcell.models import Cell, Cell_Batch
    from applib.plate.multimode_reader import multimodereader_xls
    from applib.data.set_fielddata import set_Fields_fromDict
    from dscreen.models import Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")


        


   # TestPlate XLSX -------------------------------------------------------------
    if prgArgs.table == 'TestPlateList':
        if prgArgs.excelfile and prgArgs.runid:

            OutNumbers = {'Processed Assays':0,'New Assays':0, 'Uploaded Assays':0,
                          'Processed Plates':0,'New Plates':0, 'Uploaded Plates':0,
                          'Empty':0}
            PrepSheets = get_PlatePrep_xlsx(prgArgs.excelfile,sheets=['TestPlateList','Assays'])

            # Assays ------------------------------------------------------------------------
            AssayDict = {}
            Assay_FieldList = ['assay_subtype','test_media', 'test_dye', 'test_enviroment', 'test_time',
                                'test_temperature', 'subculture_type', 'test_addition', ]
            Assay_FKeyDict  = {'organism_id': Organism,'cell_id': Cell}

            validStatus = True
            for idx,row in tqdm(PrepSheets['Assays'].iterrows(), total= len(PrepSheets['Assays']), desc='[Assays]'):
                OutNumbers['Processed Assays'] += 1
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
                    validStatus = set_Fields_fromDict(djAss,row,FieldList=Assay_FieldList, fkeyDict=Assay_FKeyDict)
                    new_assay = True
                    OutNumbers['New Assays'] += 1

                if prgArgs.upload and new_assay and validStatus:
                    OutNumbers['Uploaded Assays'] += 1
                    djAss.save()

                AssayDict[row['assay_id']] = djAss

            logger.info(f"[Assays]: New Assays: {OutNumbers['New Assays']} of {OutNumbers['Processed Assays']} - Uploaded {OutNumbers['Uploaded Assays']}")

            # TestPlates  ------------------------------------------------------------------------
            TestPlateDict = {}
            TestPlate_FieldList = ['plating','test_media', 'test_dye','test_additive', 'processing', 'issues' ]
            TestPlate_FKeyDict  = {'assay_id': Assay,'cellbatch_id': Cell_Batch, 'orgbatch_id': Organism_Batch, 'labware_id': Labware}
            TestPlate_DictList = ['result_type']

            for idx,row in tqdm(PrepSheets['TestPlateList'].iterrows(), total= len(PrepSheets['TestPlateList']), desc='[TestPlates]'):
                OutNumbers['Processed Plates'] += 1
                validStatus = True
                djTP = TestPlate.get(row['plate_id'],WellData=False)
                if djTP is None:
                    logger.info(f"[TestPlate] {row['plate_id']} does not exist - Upload first the ReadOuts or check the Plate_ID")
                    OutNumbers['New Plates'] += 1
                else:
                    validStatus = set_Fields_fromDict(djTP,row,FieldList=TestPlate_FieldList, ArrayDict={}, DictList=TestPlate_DictList, fkeyDict=TestPlate_FKeyDict)
                    TestPlateDict[row['plate_id']] = djTP

                if prgArgs.upload and validStatus:
                    OutNumbers['Uploaded Assays'] += 1
                    djTP.save()

            logger.info(f"[TestPlates]: New Plates: {OutNumbers['New Plates']} of {OutNumbers['Processed Plates']} - Uploaded {OutNumbers['Uploaded Plates']}")

            if not validStatus:
                logger.info(f"[TestPlates]: New Plates: {OutNumbers['New Plates']} of {OutNumbers['Processed Plates']} - Uploaded {OutNumbers['Uploaded Plates']}")

            logger.info(f" [TestPlateList] {OutNumbers}")
    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [TestPlate]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
    prgParser.add_argument("-e","--excel",default=None,required=True, dest="excelfile", action='store', help="Excel File")
    prgParser.add_argument("--prefix",default=None,required=False, dest="prefix", action='store', help="Prefix to add to PlateID")

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