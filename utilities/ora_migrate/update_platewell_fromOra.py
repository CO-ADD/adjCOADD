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
logDir = "log"
logFileName = os.path.join(logDir,f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

if not os.path.isdir(logDir):
    os.mkdir(logDir)

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

#def convert_oraCmpBatch_djCmpBatch(CmpLst):

def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from applib.data.set_fielddata import set_arrayFields, set_Fields, set_Dictionaries, set_fkeyFields, set_arrayDictionaries
    from dplate.models import Labware, TestPlate, TestWell,MasterPlate, MasterWell
    from dsample.models import Convert_ProjectID, Convert_CompoundID
    from dscreen.models import Screen_Run
    from dorganism.models import Organism_Batch
    from dcell.models import Cell_Batch
    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID
    from update_utils import convert_castdb_compoundid_from_ora

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # For Testing
    # D000092208 (AntiBio_R009) - DR compound_id is null
    # E00095668 (OPXP_R25) - synMIC set_id/conc_type not null
    # HC162-09-21 
    # E00092682 - SC 

    # TestWells -------------------------------------------------------------
    if prgArgs.table == "TestWells" :

        OutName = "[TestWells]"
        OutDict = []
        OutFile = f"UpdateTestWells_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        logger.info(f"{OutName} ---------------------------------------------------------")
        CastDB = openCastDB()

        twSQL = ""
        if prgArgs.plateid:
            twSQL = f"Select * From TestWell Where plate_id = '{prgArgs.plateid}'"
            nWells = CastDB.nCount(f"Select count(1) From TestWell Where plate_id = '{prgArgs.plateid}'" )
        elif prgArgs.new:
            twSQL = "Select * From TestWell Where is_migrated < 1"
            nWells = CastDB.nCount("Select count(1) From TestWell Where is_migrated < 1" )

        elif int(prgArgs.test) > 0:
            twSQL += f" Fetch First {int(prgArgs.test)} Rows Only "
            nWells = int(prgArgs.test)
        else:
            twSQL = "Select * From TestWell "
            nWells = CastDB.nCount("Select count(1) From TestWell" )
        logger.info(f"{OutName} {nWells} ")

        # ---------------------------------------------------------------------------------
        # Settings
        renameCol = {
            "iscontrol":  "is_control",
            "isposcontrol":    "is_poscontrol",
            "isnegcontrol":    "is_negcontrol",
            "issample":       "is_sample",
            "isvalid":      "is_valid",
            "isskip":   "is_skip",
            "active":   "act_type",
        }

        replaceValues = {
            'result_type':{'HC10':'HC50'},
            'plate_size':{'384w':384,'96w':96,},
        }

        arrayFields = {'cmpbatch_lst':['compound_id','compound2_id','compound3_id','compound4_id'],
                       'conc_lst':['conc','conc2','conc3','conc4',],
                       'conc_unit_lst':['conc_unit','conc2_unit','conc3_unit','conc4_unit'], 
                       'conc_type_lst':['conc_type','conc2_type','conc3_type','conc4_type'], 
                       'set_lst':['set_id','set2_id','set3_id','set4_id'], 
                       'readouts':['readout','readouta','readoutb'], 
                       }
        
        copyFields = ['zscore','mscore','act_type',
                      'inhibition','active','pscore',
                      'is_skip','is_sample', 'is_negcontrol', 'is_poscontrol','is_control','is_valid',
                      'volume',
                      'solvent', 'solvent_conc',
                      ]
        dictFields = ['solvent_conc_unit']
        fkeyFields = {'plate_id':TestPlate}


        CastDB.exec(twSQL)  
        sql_columns = [i[0].lower() for i in CastDB.cursor.description]
        logger.info(sql_columns)

        for crow in tqdm(CastDB.cursor, total=nWells, desc=OutName):
#        for crow in CastDB.cursor:

            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True

            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]

            for k_old in renameCol:
                row[renameCol[k_old]] = row.pop(k_old)
            if 'solvent_conc_unit' in row:
                if row['solvent_conc_unit'] is not None:
                    row['solvent_conc_unit'] = row['solvent_conc_unit'].lower()

            debugWells = ['A18','A17']
            # if row['well_id'] in debugWells:
            #     print(" [D]---------------------------------------------------------------------------------")
            # print(f" [D] {row}")

            NewEntry = False
            djPlate = TestPlate.get(row['plate_id'])
            if djPlate:
                djWell = TestWell.get(djPlate,row['well_id'])
                if djWell is None:
                    djWell = TestWell()
                    djWell.plate_id = djPlate
                    djWell.well_id = row['well_id']
                    NewEntry = True
                    OutNumbers['New Entry'] += 1

                    #print(f" {OutName} {djWell} ")

                set_Fields(djWell,row,copyFields)
                set_arrayFields(djWell,row,arrayFields)
                set_Dictionaries(djWell,row,dictFields)

                # Fix Readout_types
                _readout = []
                if len(djWell.readouts) == 1:
                    _readout.append(djPlate.readout_type)
                elif len(djWell.readouts) > 1:
                    if djPlate.readout_type == 'OD570-600':
                        _readout.append('OD570-600')
                        _readout.append('OD570')
                        _readout.append('OD600')
                djWell.readout_types = _readout

                _new_lst = []
                for _old in djWell.cmpbatch_lst:
                    _new,_valstatus = convert_castdb_compoundid_from_ora(_old)
                    _new_lst.append(_new)
                    if not _valstatus:
                        validStatus = False

                if validStatus:
                    djWell.set_cmpbatch_id(_new_lst) 

                    validDict = djWell.check_cmpbatch_id()
                    if validDict:
                        validStatus = False    
                        #print(validDict)
                        row.update(validDict)

                    validDict = djWell.check_conc_unit_dictionary()
                    if validDict:
                        validStatus = False    
                        #print(validDict)
                        row.update(validDict)

                    djWell.init_fields()
                    validDict = djWell.validate_fields(exclude=list(arrayFields.keys()))
                    if validDict:
                        validStatus = False
                        # for k in validDict:
                        #     print('Warning',k,validDict[k],'-')
                        row.update(validDict)

                    # if djWell.well_id in debugWells:
                    #     print(f" [D99] {djWell.well_id} {djWell.cmpbatch_lst}")

                    if validStatus:
                        if prgArgs.upload:
                            if NewEntry or prgArgs.overwrite:
                                djWell.chk_migration = 0
                                OutNumbers['Upload Entries'] += 1
                                djWell.save(user=prgArgs.appuser)
                                # if djWell.well_id in debugWells:
                                #     print(f" [DSAVE] {djWell.well_id} {djWell.cmpbatch_lst}")
                    else:
                        logger.info(f" [Error] Issues with {djWell.plate_id} {djWell.well_id}")
                        OutDict.append(row)
                else:
                    row.update({'Error': ' Old Compound_ID not found'})
                    OutDict.append(row)
            else:
                logger.info(f"{OutName} Plate {row['plate_id']} not found")
                OutDict.append(row)

            #print(f"{row['plate_id']} {row['well_id']}")
        if len(OutDict) > 0:
            logger.info(f"Writing Issues: {OutFile}")
            outDF = pd.DataFrame(OutDict)
            outDF.to_excel(OutFile)
        else:
            logger.info(f"No Issues")

        logger.info(f"{OutName} {OutNumbers}")
        #print(OutDict)
        CastDB.close()

        # tpDF = get_oraTestWells(int(prgArgs.test))
        logger.info("--------------------------------------------------------------------")
#        print(f"{OutName} {tpDF.columns} ")

    # MasterWells -------------------------------------------------------------
    if prgArgs.table == "MasterWells" :

        OutName = "[MasterWells]"
        OutDict = []
        OutFile = f"UpdateMasterWells_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        logger.info(f"{OutName} ---------------------------------------------------------")
        CastDB = openCastDB()

        mwSQL = ""
        if prgArgs.plateid:
            mwSQL = f"Select * From MasterWell Where plate_id = '{prgArgs.plateid}' "
            nWells = CastDB.nCount(f"Select count(1) From MasterWell Where plate_id = '{prgArgs.plateid}' " )
        elif prgArgs.new:
            mwSQL = "Select * From MasterWell Where is_migrated < 1 "
            nWells = CastDB.nCount("Select count(1) From MasterWell Where is_migrated < 1 " )

        elif int(prgArgs.test) > 0:
            mwSQL += f" Fetch First {int(prgArgs.test)} Rows Only "
            nWells = int(prgArgs.test)
        else:
            mwSQL = "Select * From MasterWell "
            nWells = CastDB.nCount("Select count(1) From Masterwell " )
        logger.info(f"{OutName} {nWells} ")

        # ---------------------------------------------------------------------------------
        # Settings
        renameCol = {
            # "iscontrol":  "is_control",
            # "isposcontrol":    "is_poscontrol",
            # "isnegcontrol":    "is_negcontrol",
            # "issample":       "is_sample",
            # "isvalid":      "is_valid",
            # "isskip":   "is_skip",
            # "active":   "act_type",
        }

        replaceValues = {
            # 'result_type':{'HC10':'HC50'},
            # 'plate_size':{'384w':384,'96w':96,},
        }

        arrayFields = {'cmpbatch_lst':['compound_id','compound2_id','compound3_id','compound4_id'],
                       'conc_lst':['conc','conc2','conc3','conc4',],
                       'conc_unit_lst':['conc_unit','conc2_unit','conc3_unit','conc4_unit'], 
                       'conc_type_lst':['conc_type','conc2_type','conc3_type','conc4_type'], 
                       'set_lst':['set_id','set2_id','set3_id','set4_id'], 
                       'test_conc_lst':['test_conc','test_conc2','test_conc3','test_conc4',],
                       'test_conc_unit_lst':['test_conc_unit','test_conc2_unit','test_conc3_unit','test_conc4_unit'], 
                       'test_conc_type_lst':['test_conc_type','test_conc2_type','test_conc3_type','test_conc4_type'], 
                       'dilution_lst':['dilution','dilution2','dilution3','dilution4'], 
                       }
        
        copyFields = [ 'solvent', 'solvent_conc', 'amount','volume',
                    #   'zscore','mscore','act_type',
                    #   'inhibition','active','pscore',
                    #   'is_skip','is_sample', 'is_negcontrol', 'is_poscontrol','is_control','is_valid',
                    #   'volume',
                    #   'solvent', 'solvent_conc',
                      ]
        dictFields = ['solvent_conc_unit','amount_unit','volume_unit']
        fkeyFields = {'plate_id':MasterPlate}


        CastDB.exec(mwSQL)  
        sql_columns = [i[0].lower() for i in CastDB.cursor.description]
        logger.info(sql_columns)

        for crow in tqdm(CastDB.cursor, total=nWells, desc=OutName):
#        for crow in CastDB.cursor:
            #print(f"{crow}")
            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True

            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]

            for k_old in renameCol:
                row[renameCol[k_old]] = row.pop(k_old)

            debugWells = ['A18','A17']
            # if row['well_id'] in debugWells:
            #     print(" [D]---------------------------------------------------------------------------------")
            # print(f" [D] {row}")

            NewEntry = False
            djPlate = MasterPlate.get(row['plate_id'])
            if djPlate:
                djWell = MasterWell.get(djPlate,row['well_id'])
                if djWell is None:
                    djWell = MasterWell()
                    djWell.plate_id = djPlate
                    djWell.well_id = row['well_id']
                    NewEntry = True
                    OutNumbers['New Entry'] += 1

                    #print(f" {OutName} {djWell} ")

                set_Fields(djWell,row,copyFields)
                set_arrayFields(djWell,row,arrayFields)
                set_Dictionaries(djWell,row,dictFields)

                _new_lst = []
                for _old in djWell.cmpbatch_lst:
                    _new,_valstatus = convert_castdb_compoundid_from_ora(_old)
                    _new_lst.append(_new)
                    if not _valstatus:
                        validStatus = False

                if validStatus:
                    djWell.set_cmpbatch_id(_new_lst) 

                    validDict = djWell.check_cmpbatch_id()
                    if validDict:
                        validStatus = False    
                        #print(validDict)
                        row.update(validDict)

                    validDict = djWell.check_conc_unit_dictionary()
                    if validDict:
                        validStatus = False    
                        #print(validDict)
                        row.update(validDict)

                    djWell.init_fields()
                    validDict = djWell.validate_fields(exclude=list(arrayFields.keys()))
                    if validDict:
                        validStatus = False
                        # for k in validDict:
                        #     print('Warning',k,validDict[k],'-')
                        row.update(validDict)

                    # Fix Barcode '' -> Null
                    if hasattr(djWell,'barcode'):
                        if djWell.barcode == '':
                            djWell.barcode = None

                    # if djWell.well_id in debugWells:
                    #     print(f" [D99] {djWell.well_id} {djWell.cmpbatch_lst}")

                    if validStatus:
                        if prgArgs.upload:
                            if NewEntry or prgArgs.overwrite:
                                djWell.chk_migration = 0
                                OutNumbers['Upload Entries'] += 1
                                djWell.save(user=prgArgs.appuser)
                                # if djWell.well_id in debugWells:
                                #     print(f" [DSAVE] {djWell.well_id} {djWell.cmpbatch_lst}")
                    else:
                        logger.info(f" [Error] Issues with {djWell.plate_id} {djWell.well_id}")
                        OutDict.append(row)
                else:
                    row.update({'Error': ' Old Compound_ID not found'})
                    OutDict.append(row)
            else:
                logger.info(f"{OutName} Plate {row['plate_id']} not found")
                OutDict.append(row)

            #print(f"{row['plate_id']} {row['well_id']}")
        if len(OutDict) > 0:
            logger.info(f"Writing Issues: {OutFile}")
            outDF = pd.DataFrame(OutDict)
            outDF.to_excel(OutFile)
        else:
            logger.info(f"No Issues")

        logger.info(f"{OutName} {OutNumbers}")
        #print(OutDict)
        CastDB.close()

        logger.info("--------------------------------------------------------------------")
#        print(f"{OutName} {tpDF.columns} ")

        # ------------------------------


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

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)


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
