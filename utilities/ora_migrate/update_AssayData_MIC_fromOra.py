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
logName = "Upload_AssayMIC"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from adjCOADD.applib.data.set_fielddata import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields, set_arrayDictionaries
    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import Convert_ProjectID, Convert_CompoundID
    from dscreen.models import AssayData_MIC, Screen_Run
    from dorganism.models import Organism_Batch
    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID
    from update_utils import convert_castdb_compoundid_from_ora
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == "AssayMIC" :

        OutName = "[AssayData_MIC]"
        OutDict = []
        OutFile = f"UpdateAssayMIC_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        logger.info(f"{OutName} ---------------------------------------------------------")
        CastDB = openCastDB()


        assSQL = ""
        if prgArgs.plateid:
            assSQL = f"Select * From AssayData_MIC Where testplate_id = '{prgArgs.plateid}'"
            nEntries = CastDB.nCount(f"Select count(1) From AssayData_MIC Where testplate_id = '{prgArgs.plateid}'" )
        elif prgArgs.new:
            assSQL = "Select * From AssayData_MIC Where is_migrated < 1"
            nEntries = CastDB.nCount("Select count(1) From AssayData_MIC Where is_migrated < 1" )
        else:
            assSQL = "Select * From AssayData_MIC "
            nEntries = CastDB.nCount("Select count(1) From AssayData_MIC " )

        if int(prgArgs.test) > 0:
            assSQL += f" Fetch First {int(prgArgs.test)} Rows Only "
            nEntries = int(prgArgs.test)
        logger.info(f"{OutName} {nEntries} ")


        CastDB.exec(assSQL)  
        sql_columns = [i[0].lower() for i in CastDB.cursor.description]
        logger.info(sql_columns)

        # Settings -------------------------------------------------------------------
        renameCol = {
            "dmax":  "inhibit_max",
            "dmin":    "inhibit_min",
            "cmax":    "conc_max",
            "cmin":       "conc_min",
            "chk_mic":      "ref_mic_check",
            "assaytype_id":   "assay_id",
        }

        replaceValues = {
            'result_type':{'HC10':'HC50'},
            'plate_size':{'384w':384,'96w':96,},
        }
        
        arrayFields = {'cmpbatch_lst':['compound1_id','compound2_id','compound3_id','compound4_id'],
                       }
        
        copyFields = ['analysis',
                      'mic','mic_unit','mic_skips','active','act_score','pscore',
                      'data_quality','valid',
                      'inhibit_max','inhibit_min','conc_max','conc_min',
    #                  'media','dye', 'additive', 
                      'ref_mic','ref_mic_check',
                      'ic50','ic50_unit','ic50_pscore','ic50_quality','ic50_r2','ic50_slope',
                      'pub_status','pub_date',
                      ]
        
        dictFields = ['solvent_conc_unit']

        fkeyFields = {'run_id':Screen_Run}

        for crow in tqdm(CastDB.cursor, total=nEntries, desc=OutName):
#        for crow in CastDB.cursor:

            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True

            row = dict()
            for col in sql_columns:
                row[col.lower()] = crow[sql_columns.index(col)]

            for k_old in renameCol:
                row[renameCol[k_old]] = row.pop(k_old)

            NewEntry = False
            djPlate = TestPlate.get(row['testplate_id'])
            if djPlate:
                djAssay = AssayData_MIC.get(djPlate, row['testwell_id'],verbose=0)
                if djAssay is None:
                    djAssay = AssayData_MIC()
                    djAssay.testplate_id = djPlate
                    djAssay.testwell_id = row['testwell_id']
                    NewEntry = True
                    OutNumbers['New Entry'] += 1

                    #print(f" {OutName} {djWell} ")
                djAssay.assay_id = reformat_OrganismID(row['assay_id'])
                set_dictFields(djAssay,row,copyFields)
                set_arrayFields(djAssay,row,arrayFields)
                set_fkeyFields(djAssay,row,fkeyFields)
                #set_Dictionaries(djAssay,row,dictFields)

                # Fix CmpBatch ID's
                _new_lst = []
                for _old in djAssay.cmpbatch_lst:
                    _new,_valstatus = convert_castdb_compoundid_from_ora(_old)
                    _new_lst.append(_new)
                    if not _valstatus:
                        validStatus = False
                if validStatus:
                    djAssay.cmpbatch_lst = _new_lst
                    djAssay.n_cmpbatches = len(_new_lst)

                    validDict = djAssay.check_cmpbatch_id()
                    if validDict:
                        validStatus = False    
                        row.update(validDict)

                    djAssay.init_fields()
                    #validDict = djAssay.validate_fields()
                    validDict = djAssay.validate_fields(exclude=list(arrayFields.keys()))
                    if validDict:
                        validStatus = False
                        # for k in validDict:
                        #     logger.warning('Warning',k,validDict[k],'-')
                        row.update(validDict)

                    if validStatus:
                        if prgArgs.upload:
                            if NewEntry or prgArgs.overwrite:
                                djAssay.chk_migration = 0
                                OutNumbers['Upload Entries'] += 1
                                djAssay.save(user=prgArgs.appuser)
                    else:
                        logger.info(f" [Error] Issues with {djAssay.testplate_id} {djAssay.testwell_id}")
                        OutDict.append(row)
                else:
                    row.update({'Error': ' Old Compound_ID not found'})
                    OutDict.append(row)
            else:
                logger.info(f"{OutName} Plate {row['testplate_id']} not found")
                OutDict.append(row)

        logger.info(f"{OutName} {OutNumbers}")
        CastDB.close()

        if len(OutDict) > 0:
            logger.info(f"Writing Issues: {OutFile}")
            outDF = pd.DataFrame(OutDict)
            
            with pd.ExcelWriter(OutFile) as writer:
                outDF.to_excel(writer, sheet_name='Issues')


        # for e in OutDict
        #     logger.info(e)
        logger.info("--------------------------------------------------------------------")


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