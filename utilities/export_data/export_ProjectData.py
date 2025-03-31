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
logName = "CalcTestPlateDoseresponse"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
#    format="[%(name)-20s] %(message)s ",
    format="%(message)s",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)



#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    django.setup()

    from dscreen.models import Assay
    from dplate.models import Labware, TestPlate, TestWell
    from dorganism.models import Organism, Organism_Batch
    from dsample.models import Compound_Batch,COADD_Compound, ABase_Compound, Project
    from dsummary.utils.summary_data import sum_cmpbatch_sc, sum_cmpbatch_dr
    from dscreen.models import AssayData_MIC,AssayData_CC50,AssayData_HC50
    from dcell.models import Cell, Cell_Batch
    from applib.plate.multimode_reader import multimodereader_xls
    from applib.bio.doseresponse import DoseResponse
    from applib.data.set_fielddata import set_model_from_dict
    from dscreen.models import Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")


    djPrj = Project.get(prgArgs.projectid)
    if djPrj is not None:
        qryCOADD = COADD_Compound.objects.filter(project_id = prgArgs.projectid).values('compound_id','compound_code',)
#        lstCOADD = list(qryCOADD)
        dictCOADD = {}
        lstCmpBatches = []
        for qry in qryCOADD:
            dictCOADD[qry['compound_id']] = qry
            dictCOADD[qry['compound_id']]['Source'] = 'COADD'
            lstCmpBatches.append(qry['compound_id'])


        logger.info(f" [Project] {prgArgs.projectid} COADD : {len(dictCOADD)}")
        twCmp = TestWell.objects.filter(cmpbatch_lst__overlap=lstCmpBatches, 
                                        plate_id__result_type='Inhibition',
                                        plate_id__plate_quality='Valid').order_by('cmpbatch_lst').distinct('cmpbatch_lst').values_list('cmpbatch_lst')
        
        logger.info(f" [Sum CmpBatch SC] TestWell: {twCmp.count()}  ")
        for cmps in tqdm(twCmp, desc='[Sum SC]'):
            _numbers,_outdict  = sum_cmpbatch_sc(cmps[0],upload=True,overwrite=True,appuser=prgArgs.appuser)
        for cmps in tqdm(twCmp, desc='[Sum DR]'):
            _numbers,_outdict  = sum_cmpbatch_dr(cmps[0],upload=True,overwrite=True,appuser=prgArgs.appuser)



    else:
        logger.info(f" [Project] {prgArgs.projectid} NOT FOUND")




#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    # prgParser.add_argument("-t",default=None,required=False, dest="table", action='store', help="Table to upload [TestPlate]")
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    # prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-p","--project",default=None,required=True, dest="projectid", action='store', help="Project")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="RunID")
#    prgParser.add_argument("-e","--excel",default=None,required=True, dest="excelfile", action='store', help="Excel File")
    prgParser.add_argument("--plotdir",default=None,required=False, dest="plotdir", action='store', help="Folder for Plots")
    #prgParser.add_argument("-o","--outdir",default=None,required=False, dest="outdir", action='store', help="Prefix to add to PlateID")

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