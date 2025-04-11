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
                PlatePrep_Sheets[key] = pd.read_excel(xls, key).fillna('-')
                PlatePrep_Sheets[key].columns = [c.lower() for c in PlatePrep_Sheets[key].columns]
                
    return(PlatePrep_Sheets)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    django.setup()

    from dscreen.models import Assay
    from dplate.models import Labware, TestPlate, TestWell, MasterPlate
    from dorganism.models import Organism, Organism_Batch
    from dcell.models import Cell, Cell_Batch
    from applib.plate.multimode_reader import multimodereader_xls
    from applib.data.set_fielddata import set_model_from_dict
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
            ass_Fields = ['assay_subtype','sum_assay_id','test_media', 'test_dye', 'test_enviroment', 'test_time',
                                'test_temperature', 'subculture_type', 'test_addition', ]
            ass_FKeys  = {'organism_id': Organism,'cell_id': Cell}

            validStatus = True
            
            logger.info(f"[Assays]")
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
                    validStatus = set_model_from_dict(djAss,row,list_Fields=ass_Fields, dict_FKeys=ass_FKeys)
                    new_assay = True
                    OutNumbers['New Assays'] += 1

                if prgArgs.upload and new_assay and validStatus:
                    OutNumbers['Uploaded Assays'] += 1
                    djAss.save()

                AssayDict[row['assay_id']] = djAss
            
            logger.info(f"[Assays]: {OutNumbers['New Assays']} new assays (of {OutNumbers['Processed Assays']}) ")
            logger.info(f"[Assays]: New assays uploaded {OutNumbers['Uploaded Assays']} [Upload: {prgArgs.upload}]")
            if (OutNumbers['Processed Assays'] - OutNumbers['New Assays']) > 0:
                logger.info(f"[Assays]: -- Existing assays need to updated online")
            logger.info(f"[Assays]")

            # TestPlates  ------------------------------------------------------------------------
            PrepSheets['TestPlateList'].rename(columns={"test_strain": "test_orgbatch_id", "test_cell": "test_cellbatch_id"},inplace=True)
            print( PrepSheets['TestPlateList'].columns)
            TestPlateDict = {}
            tp_Fields = ['plating','test_media', 'test_dye','test_additive', 'processing', 'issues','control_layout' ]
            tp_FKeys  = {'assay_id': Assay,'test_cellbatch_id': Cell_Batch, 'test_orgbatch_id': Organism_Batch, 'labware_id': Labware}
            tp_Dicts = ['result_type']
            tp_Arrays = {'motherplate_ids':['motherplate_id','motherplate2_id'],'synergy_cmpbatches':['syn_compounds_ab','syn_compounds_pot']}
            tp_FKeyArrays = {'motherplate_ids':{'model':MasterPlate, 'fields': ['motherplate_id','motherplate2_id']},
                             'synergy_cmpbatches':{'model':Organism_Batch, 'fields':['syn_compounds_ab','syn_compounds_pot']}
                            }
            #
            # TODO check that essential columns are in Excel Sheet
            #   ['testplate_id', 'assay_id', 'test_strain','result_type','labware_id', 'control_layout' ] 
            #

            # SynCompounds
            for idx,row in tqdm(PrepSheets['TestPlateList'].iterrows(), total= len(PrepSheets['TestPlateList']), desc='[TestPlates]'):
                OutNumbers['Processed Plates'] += 1
                validStatus = True
                djTP = TestPlate.get(row['testplate_id'],WellData=False)
                if djTP is None:
                    logger.info(f"[TestPlate] {row['testplate_id']} does not exist - Upload first the ReadOuts or check the Plate_ID")
                    OutNumbers['New Plates'] += 1
                else:
                    validStatus = set_model_from_dict(djTP,row,
                                                      list_Fields=tp_Fields, 
                                                      dict_Arrays=tp_Arrays, 
                                                      list_Dicts=tp_Dicts, 
                                                      dict_FKeys=tp_FKeys,
                                                      dict_FKeyArrays= tp_FKeyArrays)
                    TestPlateDict[row['testplate_id']] = djTP
                    #print(f" [{djTP.plate_id}] {validStatus}")
                if prgArgs.upload and validStatus:
                    OutNumbers['Uploaded Plates'] += 1
                    djTP.test_orgbatch_id
                    djTP.save()

            logger.info(f"[TestPlates]")
            logger.info(f"[TestPlates]: {OutNumbers['Processed Plates']} TestPlates ({OutNumbers['New Plates']} new plates)")
            logger.info(f"[TestPlates]: Plates uploaded {OutNumbers['Uploaded Plates']} [Upload: {prgArgs.upload}]")
            if OutNumbers['New Plates'] > 0:     
                logger.info(f"[TestPlates]: -- The {OutNumbers['New Plates']} new plates are not uploaded")
                logger.info(f"[TestPlates]: -- ReadOuts need to be uploaded first, or check the TestPlate_ID")
            logger.info(f"[TestPlates]")

            # if not validStatus:
            #     logger.info(f"[TestPlates]: New Plates: {OutNumbers['New Plates']} of {OutNumbers['Processed Plates']} - Uploaded {OutNumbers['Uploaded Plates']}")

            # logger.info(f" [TestPlateList] {OutNumbers}")
    
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