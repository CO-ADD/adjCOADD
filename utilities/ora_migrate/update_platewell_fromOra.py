#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
from oraCastDB.oraCastDB import openCastDB

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_Well"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields
    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import Convert_ProjectID, Convert_CompoundID
    from dscreen.models import Screen_Run
    from dorganism.models import Organism_Batch
    from dcell.models import Cell_Batch
    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

 
   # TestWells -------------------------------------------------------------
    if prgArgs.table == "TestWells" :

        OutName = "[TestWells]"
        OutDict = []
        OutFile = f"UpdateTestWells_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        print(f"{OutName} ---------------------------------------------------------")
        CastDB = openCastDB()

        twSQL = "Select * From TestWell "
        if int(prgArgs.test) > 0:
            twSQL += f" Fetch First {int(prgArgs.test)} Rows Only "
            nWells = int(prgArgs.test)
        else:
            nWells = CastDB.nCount("Select count(1) From TestWell" )
        print(f"{OutName} {nWells} ")

        # Settings
        renameCol = {
            "media_id":  "test_media",
            "issues":    "test_issues",
            "volume":    "test_volume",
            "processing":       "test_processing",
            "n_dr":      "n_doseresponse",
            "n_syn":   "n_synergies",
            "has_readout":   "n_reads",
            "has_compound":   "n_sample",
            "has_layout":   "n_layout",
            "nreads":   "n_readouts",
            "readout_id" : "readout_type",
            "layout_control" : "control_layout",
            "assaytype_id" : "assay_id" 
        }

        replaceValues = {
            'result_type':{'HC10':'HC50'},
            'plate_size':{'384w':384,'96w':96,},
        }

        arrayFields = {'motherplate_ids':['motherplate_id','motherplate2_id'],
                       'synergy_cmpbatches':['syn_compounds_a', 'syn_compounds_b'],
                       'poscontrol_stats':['poscontrol_median','poscontrol_mad','poscontrol_ave','poscontrol_stdev'], 
                       'negcontrol_stats':['negcontrol_median','negcontrol_mad','negcontrol_ave','negcontrol_stdev'], 
                       'sample_stats':['sample_median','sample_mad','sample_ave','sample_stdev'], 
                       'edge_stats':['edge_median','nonedge_median'], 
                       }
        copyFields = ['plating','process_status',
                      'test_date','test_media','test_dye','test_additive','test_issues','test_volume',
                      'reader','experiment', 'protocol', 'inputfile','test_processing',
                      'control_layout','readout_type',
                      'plate_qc', 'zfactor','analysis_parameter',
                      'n_inhibition','n_doseresponse','n_synergies','n_readouts',
                      'n_reads', 'n_sample', 'n_layout',
                      ]
        dictFields = ['result_type','plate_quality','plate_type']
        fkeyFields = {'plate_id':TestPlate}


        CastDB.exec(twSQL)  
        sql_columns = [i[0] for i in CastDB.cursor.description]

        for crow in tqdm(CastDB.cursor, total=nWells, desc=OutName):
            # gen Dict
            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]

            # for k_old in renameCol:
            #     row[renameCol[k_old]] = row.pop(k_old)

            #print(f"{row['plate_id']} {row['well_id']}")
        CastDB.close()

        # tpDF = get_oraTestWells(int(prgArgs.test))
        print("--------------------------------------------------------------------")
#        print(f"{OutName} {tpDF.columns} ")


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

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
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
