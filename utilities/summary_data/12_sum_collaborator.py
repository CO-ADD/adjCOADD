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
logName = "SumProject"
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
    if prgArgs.collaborator:
        cAnalysis.qry_by_Collaborator(prgArgs.collaborator)
        #cAnalysis.qry_by_ProjectID(prgArgs.projectid)
        if cAnalysis.n_compounds>0:
            cAnalysis.get_dataframe(SC_Only=prgArgs.sc_only, DR_Only=prgArgs.dr_only)
            cAnalysis.get_sample_info()
            cAnalysis.get_assay_info()
            cAnalysis.get_testplate_info()
            cAnalysis.gen_pivot_tables(PivColumns=['assay_org','assay_type','result_type'])
            # # cAnalysis.add_Vitek_AST()

            cAnalysis.to_excel(prgArgs.excelfile)


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
    prgParser.add_argument("-c","--collab",default=None,required=True, dest="collaborator", action='store', help="Collaborator Code or ID")
    prgParser.add_argument("-f","--format",default=None,required=False, dest="format", action='store', help="Report Format [COADD/Total]")
#    prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
    prgParser.add_argument("-e","--excel",default=None,required=False, dest="excelfile", action='store', help="Excel File")
    prgParser.add_argument("--plotdir",default=None,required=False, dest="plotdir", action='store', help="Folder for Plots")
    #prgParser.add_argument("-o","--outdir",default=None,required=False, dest="outdir", action='store', help="Prefix to add to PlateID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    #prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    prgParser.add_argument("--sc",default=False,required=False, dest="sc_only", action='store_true', help="Only Single Concentration data")
    prgParser.add_argument("--dr",default=False,required=False, dest="dr_only", action='store_true', help="Only DoseResponse data")

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