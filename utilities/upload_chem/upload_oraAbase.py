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

# from zUtils import zData
from oraABase.oraABase import openABase, get_CompoundBatch
#from rdkit import Chem 
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ABase"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def get_AbaseRegView(test=0):

    renameCol = {
        "labware_addon":  "labware_notes",
        "material":       "plate_material",
        "work_volume" :   "working_volume",
    }

    replaceValues = {
      'plate_size':{'384w':384,'96w':96,},
    }

    _SQL = "Select * from ObjdRgst_View "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        _SQL += f" Fetch First {test} Rows Only "

    ABaseDB = openABase()
    logger.info(f"[RegView] ... ")
    _DF = pd.DataFrame(ABaseDB.get_dict_list(_SQL))
    nTotal = len(_DF)
    logger.info(f"[RegView] {nTotal} ")
    ABaseDB.close()

    # logger.info(f"DF - Rename Columns {len(renameCol)}")
    # _DF.rename(columns=renameCol, inplace=True)

    # logger.info(f"DF - Replace Values {len(replaceValues)}")
    # for k in replaceValues:
    #     _DF[k].replace(replaceValues[k],inplace=True)

    return(_DF)


def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from adjCOADD.applib.data.set_fielddata import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields
    from dsample.models import Compound_Batch
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

 
   # ABase -------------------------------------------------------------
    if prgArgs.table == "RegView" :

        OutName = "[RegView]"
        OutDict = []
        OutFile = f"UpdateRegView_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        print(f"{OutName} ---------------------------------------------------------")
        regDF = get_AbaseRegView(int(prgArgs.test))
        print("--------------------------------------------------------------------")
        print(f"{OutName} {regDF.columns} ")

        
        for idx,row in tqdm(regDF.iterrows(), total=regDF.shape[0], desc=OutName):
            #print(row)
            OutNumbers['Processed'] += 1
            NewEntry = False
            validStatus = True
            oraBatch_id = f"{row['objdid']}:{row['objdbatchref']}"
            djBatch_id = f"{row['objdid'].replace('_','')}_{row['objdbatchref']}"


            djObj = Compound_Batch.get(djBatch_id)
            if djObj is None:
                NewEntry = True
                djObj = Compound_Batch()
                djObj.cmpbatch_id = djBatch_id
                djObj.batch_id = row['objdbatchref']

            djObj.full_mf = row['rgstfullmolformula']
            djObj.full_mw = row['rgstfullmolmassvalue']
            djObj.batch_source = 'ABASE'
            djObj.batch_code = f"{row['objdid']}:{row['objdbatchref']}"
            if 'rgstdrugname' in row:
                djObj.batch_notes = row['rgstdrugname']

            djObj.init_fields()
            validDict = djObj.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('Warning',k,validDict[k],'-')
                OutDict.append(row)
            #print(f" {validStatus} {prgArgs.upload}")
            if validStatus:
                if prgArgs.upload:
                    if NewEntry or prgArgs.overwrite:
                        OutNumbers['Upload Entries'] += 1
                        djObj.save(user=prgArgs.appuser)

        print(f"{OutName} {OutNumbers}")
        print(OutDict)

 


        #get_CompoundBatch 

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