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
logName = "Update_Antibiogram"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
def main(prgArgs,djDir):
#-----------------------------------------------------------------------------------

    django.setup()

    from apputil.models import Dictionary
    from dsample.models import COADD_Compound, Compound_Batch
    from ddrug.models import Drug, MIC_COADD
    from dsummary.utils.summary_data import sum_structure_sc
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from adjcoadd.constants import COMPOUND_SEP
    from django.db.models import Q

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'MIC_COADD':
        if prgArgs.runid:
            qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                                  testplate_id__plate_quality = 'Valid',
                                                  run_id=prgArgs.runid,
                                                  )
        else:
            qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                                  testplate_id__plate_quality = 'Valid'
                                                  )
        nMIC = qryMIC.count()
        logger.info(f" [{prgArgs.table}] {nMIC} for {prgArgs.runid}")
        if qryMIC:

            OutNumbers = {'Processed':0,'New':0, 'Uploaded':0,'Empty':0, 'Failed':0}

            for mic in tqdm(qryMIC, total=nMIC, desc=prgArgs.table):
                validStatus = True 
                OutNumbers['Processed'] += 1
                cmp_lst=[]
                for batch in mic.cmpbatch_lst:
                    _l = batch.split('_')
                    cmp_lst.append("_".join(_l[:2]))
                cmps = "|".join(cmp_lst)

                if Drug.objects.filter(uq_imb=cmps).exists():
                    djDrug =   Drug.objects.get(uq_imb=cmps)               
                    djMIC = MIC_COADD.get(mic.testplate_id.test_orgbatch_id,djDrug,mic.testplate_id,mic.testwell_id,verbose=0)
                    if djMIC is None:
                        OutNumbers['New'] += 1

                        djMIC = MIC_COADD()
                        djMIC.orgbatch_id = mic.testplate_id.test_orgbatch_id
                        djMIC.drug_id = djDrug
                        djMIC.testplate_id = mic.testplate_id
                        djMIC.testwell_id = mic.testwell_id
                    djMIC.run_id = mic.run_id

                    djMIC.mic = mic.mic
                    djMIC.mic_unit = mic.mic_unit
                    djMIC.mic_type = Dictionary.get(MIC_COADD.DICTIONARY_FIELDS["mic_type"],'BMD',None,verbose=1)

                    djMIC.plate_size = Dictionary.get(MIC_COADD.DICTIONARY_FIELDS["plate_size"],mic.testplate_id.labware_id.plate_size,None,verbose=1) 
                    djMIC.plate_material = mic.testplate_id.labware_id.plate_material 

                    djMIC.set_defaults_model()
                    validDict = djMIC.validate_fields()
                    if validDict:
                        validStatus = False
                        OutNumbers['Failed'] += 1 
                        logger.error(f"[{prgArgs.table}] {validDict}")

                    if prgArgs.upload and validStatus:
                        djMIC.save()
                        OutNumbers['Uploaded'] += 1
                else:
                    logger.warning(f" [MIC COADD] Drug {cmps} not found")

        logger.info(f" [{prgArgs.table}] {OutNumbers}")





#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [MIC_COADD]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

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