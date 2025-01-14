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
logName = "Upload_Assay_Xlsx"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

        
#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from dsummary.utils.upd_sum_cmpbatch import sum_cmpbatch_sc
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run, Assay
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from apputil.utils.set_data import set_dictFields
    from dcell.models import Cell
    from dorganism.models import Organism
    from adjcoadd.constants import COMPOUND_SEP
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # AssayData MIC -------------------------------------------------------------

    choiceTables = ['Assay']
    if prgArgs.table in choiceTables:

        ExcelFile = "C:/Data/CastDB/Migration/AssayData/TestPlate_AssayID_v01.xlsx"
        SheetName = "Assays"
        OutFile = "UploadProject_Issues.xlsx"

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

        if prgArgs.table == 'Assay':
            print(f"Reading {ExcelFile}.[{SheetName}]")
            df = pd.read_excel(ExcelFile,sheet_name=SheetName).fillna('')
            print(df.columns)

            OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0,'Empty Entries':0}
            OutDict = []
            # # CONTACT_A_ID	CONTACT_B_ID	
            # ignoreFields = ['COUNTRY','STATUS','PROJECT_ACTION','SCREEN_CONC','SCREEN_CONC_UNIT','COADD_ID','ANTIMICRO_STATUS']
            cpyFields = ['ora_assay_id','sum_assay_id',
                        'assay_panel','assay_code','assay_type','assay_subtype',
                        'test_media','test_dye','test_additive'
                        ]
            # arrayFields = {'screen_status': 'screen_status',
            #                'report_status': 'report_status',
            #                'compound_status': 'compound_status',
            #                'data_status': 'data_status',
            #                'stock_status': 'stock_status',
            #                'pub_status': 'pub_status',
            #                'ora_contact_ids':['CONTACT_A_ID','CONTACT_B_ID']
            #             }
            # dictFields = ['project_type','provided_container','stock_conc_unit',]
                
            for idx,row in tqdm(df.iterrows(), total=df.shape[0]):
                NewEntry = False
                validStatus = True

                djAss = Assay.get(row['assay_id'])
                if djAss is None:
                    NewEntry = True
                    djAss = Assay()
                    djAss.assay_id = row['assay_id']
                        
                set_dictFields(djAss,row,cpyFields)

                if row['organism_id'] != '' :
                    djAss.organism_id = Organism.get(row['organism_id'])
                # else:
                #     djAss.organism_id = None

                if row['cell_id'] != '' :
                    djAss.cell_id = Cell.get(row['cell_id'])
                # else:
                #     djAss.cell_id = None

                # Validate and Save
                djAss.init_fields()
                validDict = djAss.validate_fields()
                if validDict:
                    validStatus = False
                    # for k in validDict:
                    #     print('Warning',k,validDict[k],'-')
                    row.update(validDict)
                    logger.warning(f"{djAss.assay_id} {validDict} ")
                    OutDict.append(row)

                if validStatus:
                    if prgArgs.upload:
                        if NewEntry or prgArgs.overwrite:
                            #djSum.chk_migration = 0
                            OutNumbers['Upload Entries'] += 1
                            djAss.save(user=prgArgs.appuser)

            # if len(outDict) > 0:
            #     print(f"Writing Issues: {OutFile}")
            #     outDF = pd.DataFrame(outDict)
            #     outDF.to_excel(OutFile)
            # else:
            print(f"No Issues")


                # djPrj = Project.get(row['project_id'])
                # if djPrj is None:
                #     djPrj = Project()
                #     djPrj.project_id = row['project_id']
                #     djPrj.old_project_id = row['old_project_id']
                # if prgArgs.upload:
                #     djPrj.save()    

    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [CompoundID]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.django == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.django == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.django == 'Laptop':
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================