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
from oraCastDB import oraCastDB

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_Plate"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------
def get_oraLabware(test=0):
    from oraCastDB.oraCastDB import openCastDB

    renameCol = {
        "labware_addon":  "labware_notes",
        "material":       "plate_material",
        "work_volume" :   "working_volume",
    }

    replaceValues = {
      'plate_size':{'384w':384,'96w':96,},
    }

    lwSQL = "Select * From Labware "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        lwSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[Labware] ... ")
    lwDF = pd.DataFrame(CastDB.get_dict_list(lwSQL))
    nTotal = len(lwDF)
    logger.info(f"[Labware] {nTotal} ")
    CastDB.close()

    logger.info(f"DF - Rename Columns {len(renameCol)}")
    lwDF.rename(columns=renameCol, inplace=True)

    logger.info(f"DF - Replace Values {len(replaceValues)}")
    for k in replaceValues:
        lwDF[k].replace(replaceValues[k],inplace=True)

    return(lwDF)

#-----------------------------------------------------------------------------
def get_oraTestPlates(test=0):
    from oraCastDB.oraCastDB import openCastDB

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

    tpSQL = "Select * From TestPlate "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        tpSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[TestPlates] ... ")
    tpDF = pd.DataFrame(CastDB.get_dict_list(tpSQL))
    nTotal = len(tpDF)
    logger.info(f"[TestPlates] {nTotal} ")
    CastDB.close()

    logger.info(f"DF - Rename Columns {len(renameCol)}")
    tpDF.rename(columns=renameCol, inplace=True)

    logger.info(f"DF - Replace Values {len(replaceValues)}")
    for k in replaceValues:
        tpDF[k].replace(replaceValues[k],inplace=True)

    return(tpDF)

#-----------------------------------------------------------------------------
def get_oraTestWells(test=0):
    from oraCastDB.oraCastDB import openCastDB

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
        "assaytype_id" : "assay_id",
    }

    replaceValues = {
      'result_type':{'HC10':'HC50'},
      'plate_size':{'384w':384,'96w':96,},
    }

    twSQL = "Select * From TestWell "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        twSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[TestWells] ... ")
    twDF = pd.DataFrame(CastDB.get_dict_list(twSQL))
    nTotal = len(twDF)
    logger.info(f"[TestWells] {nTotal} ")
    CastDB.close()

    # logger.info(f"DF - Rename Columns {len(renameCol)}")
    # twDF.rename(columns=renameCol, inplace=True)

    # logger.info(f"DF - Replace Values {len(replaceValues)}")
    # for k in replaceValues:
    #     twDF[k].replace(replaceValues[k],inplace=True)

    return(twDF)


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from adjCOADD.applib.data.set_fielddata import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields
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

   # TestPlate -------------------------------------------------------------
    if prgArgs.table == "TestPlates" :

        OutName = "[TestPlates]"
        OutDict = []
        OutFile = f"UpdateTestPlates_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        logger.info(f"{OutName} ---------------------------------------------------------")
        tpDF = get_oraTestPlates(int(prgArgs.test))
        logger.info("--------------------------------------------------------------------")
        logger.info(f"{OutName} {tpDF.columns} ")

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
        fkeyFields = {'labware_id':Labware, 'run_id':Screen_Run, 'test_orgbatch':Organism_Batch, 'test_cellbatch':Cell_Batch}

        CL_replaceAssayID  = {'MA_007':['CL_0031','CL_0031_03'],   
                              'MA_008':['CL_0037','CL_0037_04'],   
                              'MA_014':['CL_0038','CL_0038_01'],   
                              'MA_021':['CL_0040','CL_0040_01'],
                              'HA_150':['CL_0078','CL_0078_01'],
                    }
        No_replaceAssayID = ['QC_LCMS','CMC_01']


        for idx,row in tqdm(tpDF.iterrows(), total=tpDF.shape[0], desc=OutName):
            #print(row)
            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True

            row['plate_type'] = 'Test'

            djObj = TestPlate.get(row['plate_id'])
            if djObj is None:
                djObj = TestPlate()
                djObj.plate_id = row['plate_id']
                
                NewEntry = True
                OutNumbers['New Entry'] += 1

            djObj.set_platesize(row['plate_size'])
            if row['assay_id'] is not None:
                if 'MA_' in row['assay_id'] or 'HA' in row['assay_id']:
                    djObj.assay_id = CL_replaceAssayID[row['assay_id']][0]
                    row['test_orgbatch'] = None
                    row['test_cellbatch'] = CL_replaceAssayID[row['assay_id']][1]
                elif row['assay_id'] in No_replaceAssayID :
                    djObj.assay_id = row['assay_id']
                    row['test_orgbatch'] = None
                    row['test_cellbatch'] = None
                else:
                    djObj.assay_id = reformat_OrganismID(row['assay_id'])
                    row['test_orgbatch'] = reformat_OrgBatchID(row['test_strain'])
                    row['test_cellbatch'] = None
            else:
                djObj.assay_id = ''
                row['test_orgbatch'] = None
                row['test_cellbatch'] = None

            set_dictFields(djObj,row,copyFields)
            set_arrayFields(djObj,row,arrayFields)
            set_Dictionaries(djObj,row,dictFields)
            set_fkeyFields(djObj,row,fkeyFields)

            djObj.init_fields()
            validDict = djObj.validate_fields(exclude=list(arrayFields.keys()))
            if validDict:
                validStatus = False
                for k in validDict:
                    logger.warning('Warning',k,validDict[k],'-')
                OutDict.append(row)

            if validStatus:
                if prgArgs.upload:
                    if NewEntry or prgArgs.overwrite:
                        OutNumbers['Upload Entries'] += 1
                        djObj.save(user=prgArgs.appuser)
        logger.info(f"{OutName} {OutNumbers}")
        logger.info(OutDict)


   # Labware -------------------------------------------------------------
    elif prgArgs.table == "Labware" :

        OutName = "[Labware]"
        OutDict = []
        OutFile = f"UpdateLabware_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        logger.info(f"{OutName} ---------------------------------------------------------")
        lwDF = get_oraLabware(int(prgArgs.test))
        logger.info("--------------------------------------------------------------------")
        logger.info(f"{OutName} {lwDF.columns} ")


        arrayFields = {}
        copyFields = ['labware_name', 'labware_type', 'labware_notes',
                'plate_size', 'brand', 'model', 
                'well_type', 'well_size', 'well_shape', 'well_bottom', 
                'working_volume','plate_color']
        dictFields = ['plate_material']

        for idx,row in tqdm(lwDF.iterrows(), total=lwDF.shape[0], desc=OutName):
            #print(row)
            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True

            djObj = Labware.get(row['labware_id'])
            if djObj is None:
                djObj = Labware()
                djObj.labware_id = row['labware_id']
                NewEntry = True
                OutNumbers['New Entry'] += 1

            set_dictFields(djObj,row,copyFields)
            set_arrayFields(djObj,row,arrayFields)
            set_Dictionaries(djObj,row,dictFields)

            djObj.init_fields()
            validDict = djObj.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    logger.warning('Warning',k,validDict[k],'-')
                OutDict.append(row)

            if validStatus:
                if prgArgs.upload:
                    if NewEntry or prgArgs.overwrite:
                        OutNumbers['Upload Entries'] += 1
                        djObj.save(user=prgArgs.appuser)
        logger.info(f"{OutName} {OutNumbers}")
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
