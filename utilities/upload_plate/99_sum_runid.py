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
logName = "SumRunID"
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
    from dsample.models import Compound_Batch
    from dcell.models import Cell, Cell_Batch
    from applib.plate.multimode_reader import multimodereader_xls
    from applib.bio.doseresponse import DoseResponse
    from applib.data.set_fielddata import set_model_from_dict
    from dsummary.utils.analyse_data import Analysis_Screening
    from dscreen.models import Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")


    cAnalysis = Analysis_Screening()
    # Process TestPlate -----------------------------------------------------------
    if prgArgs.runid:
        cAnalysis.qry_by_RunID(prgArgs.runid)
        cAnalysis.get_dataframe()
        cAnalysis.get_sample_info()
        cAnalysis.get_assay_info()

        if 'TestPlates' in prgArgs.adddata:
            cAnalysis.get_testplate_info()
        if 'Vitek' in prgArgs.adddata:
            cAnalysis.add_vitek_ast()
        if 'COADD' in prgArgs.adddata:
            cAnalysis.add_antibiogram_data(cAnalysis.ORGANISMS['COADD'])

        cAnalysis.gen_pivot_tables(PivTables = ['Values','AssayID'])

        cAnalysis.to_excel(prgArgs.excelfile, Transpose_PivTables=prgArgs.transpose)


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    # prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [TestPlate]")
    # prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    # prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    # prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    #prgParser.add_argument("-p","--pivot",default=None,required=False, dest="pivot", action='store', help="Pivot tables")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=[],required=True, dest="runid", action='store', help="RunID's as list [,]",
                type=lambda s: [item for item in s.split(',')])
    prgParser.add_argument("-a","--adddata",default=[],required=False, dest="adddata", action='store', help="Add data [Vitek/COADD] as list [,]",
                type=lambda s: [item for item in s.split(',')])
    prgParser.add_argument("-e","--excel",default=None,required=False, dest="excelfile", action='store', help="Excel File")
    # prgParser.add_argument("-f","--format",default='Check',required=False, dest="pivot", action='store', help="Format of output EXcel")
    prgParser.add_argument("--plotdir",default=None,required=False, dest="plotdir", action='store', help="Folder for Plots")
    #prgParser.add_argument("-o","--outdir",default=None,required=False, dest="outdir", action='store', help="Prefix to add to PlateID")

    prgParser.add_argument("-t",default=False,required=False, dest="transpose", action='store_true', help="Transpose pivTable")
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