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
logName = "CalcTestPlateInhibition"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
#    format="[%(name)-20s] %(message)s ",
    format="%(message)s",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    django.setup()

    from dscreen.models import Assay
    from dplate.models import Labware, TestPlate, TestWell
    from dorganism.models import Organism, Organism_Batch
    from dcell.models import Cell, Cell_Batch
    from applib.plate.multimode_reader import multimodereader_xls
    from applib.data.set_fielddata import set_model_from_dict
    from dscreen.models import Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir['djPrj']}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # TestPlate XLSX -------------------------------------------------------------
    if prgArgs.table == 'TestPlateInhibition':
        if prgArgs.runid:

            OutNumbers = {'Processed Plates':0,'Valid Plates':0, 'Rejected Plates':0, 'Failed Plates':0}

            qryTP = TestPlate.objects.filter(run_id = prgArgs.runid).values('plate_id')
            # qryTP = TestPlate.objects.filter(run_id = prgArgs.runid)
            nCnt = qryTP.count()
            logger.info(f" [{prgArgs.table}] {prgArgs.runid} : {nCnt}")

            for tp in tqdm(qryTP, desc='Testplates'):
                OutNumbers['Processed Plates'] += 1

                djTP = TestPlate.get(tp['plate_id'],WellData=True)
                if djTP.n_wells > 0 and djTP.n_reads > 0 and djTP.control_layout:
                    #print(f"{djTP.plate_id} {djTP.n_wells} {djTP.control_layout}")

                    if djTP.apply_layout() > 0:
                        djTP.calc_inhibition()
                        #print(f" [{djTP.plate_quality}]")

                        if str(djTP.plate_quality) == 'Valid':
                            OutNumbers['Valid Plates'] += 1
                        elif str(djTP.plate_quality) == 'Rejected':
                            OutNumbers['Rejected Plates'] += 1
                        else:
                            OutNumbers['Failed Plates'] += 1

                        if prgArgs.plotdir:
                            djTP.plot_heatmap('readout_1',prgArgs.plotdir)

                        if prgArgs.upload:
                            djTP.save()

                else:
                    OutNumbers['Failed Plates'] += 1
                    logger.warning(f" FAILED: {djTP.plate_id} only [Wells: {djTP.n_wells} Reads: {djTP.n_reads} or {djTP.control_layout}]")


            logger.info(f"[TestPlates]: {OutNumbers['Valid Plates']} Valid,   {OutNumbers['Rejected Plates']} Rejected, {OutNumbers['Failed Plates']} Failed of {OutNumbers['Processed Plates']} Plates")

    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [TestPlateInhibition]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")
#    prgParser.add_argument("--new",default=False,required=False, dest="new", action='store_true', help="Not migrated entries only")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--plate",default=None,required=False, dest="plateid", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=True, dest="runid", action='store', help="RunID")
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