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

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Sum_CmpBatch_Inhib"
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
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields, set_arrayDictionaries
    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from dorganism.models import Organism_Batch
#    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID
#    from update_utils import convert_castdb_compoundid_from_ora
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    
   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'Summary_CmpBatch_Inhib':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"Update{prgArgs.table}_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        qrySources = ['COADD']
        if int(prgArgs.test) > 0:
            qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')[:int(prgArgs.test)]
        else:
            qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')
        
        for cmp in qryCmp:
            CmpBatchID = cmp[0]
            qryCmpBatchID = [CmpBatchID]
            qryNCmpBatches = len(qryCmpBatchID)

            qryTW = TestWell.objects.filter(cmpbatch_lst__contains = [CmpBatchID], 
                                            n_cmpbatches = qryNCmpBatches, 
                                            plate_id__result_type = 'Inhibition'
                                           ).values(
                                               'plate_id','well_id','plate_id__result_type','plate_id__assay_id',
                                               'inhibition','mscore','act_type'
                                                 )
            DF = pd.DataFrame(qryTW)

            print(DF)                                  
            # for tw in qryTW:
            #     print(CmpBatchID, tw.plate_id, tw.well_id, tw.plate_id.result_type, tw.plate_id.assay_id, tw.inhibition)

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