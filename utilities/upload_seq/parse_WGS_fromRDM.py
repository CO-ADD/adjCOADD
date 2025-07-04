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
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_Assembly"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------

def main(prgArgs,djDir):

    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    #from applib.data.set_fielddata import set_model_arrayfields, set_dictFields, set_model_dicts
    from applib.data.str_lists import listFolders
    from apputil.utils.validation_log import Validation_Log
    
    from dgene.models import Gene,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM
    from dgene.utils.upload_gene import (WGS_RDM, get_RDM, split_BatchID_RunID, get_subdir,
                                        upload_Trim, upload_CheckM, upload_FastA, upload_AMR)
    from dgene.utils.parse_wgs import (get_FastQC_Info, get_CheckM_Info, 
                                   get_Kraken_Info, get_MLST_Info, get_GTDBTK_Info, 
                                   get_AMRFinder_Info, get_Abricate_Info, get_RGI_Info)
 
    #from dgene.utils.import_gene import (imp_Sequence_fromDict)
    #from dgene.utils.parse_wgs import ()
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

    # Assembly -------------------------------------------------------------
    nProc = {}
    nProc['Processed'] = 0
    nProc['Assembly'] = 0
    nProc['FastA'] = 0

    if (prgArgs.orgbatchid and prgArgs.runid):
        print(f"=[PARSE] Single - [{prgArgs.orgbatchid} {prgArgs.runid}]  {prgArgs.process}")
        # RDM = get_RDM(djDir['rdmDir_WGS'])
        # FastABase = os.path.join(RDM['base'],RDM['fasta'])

        WGS = WGS_RDM(djDir['rdmDir_WGS'],prgArgs.orgbatchid,prgArgs.runid,
                      SeqMethod='Illumina',
                      valLog=Validation_Log('WGS-Upload'))
        lRGI = get_RGI_Info(WGS.fasta_dir,WGS.orgbatch_id, WGS.run_id)
        nProc['Processed'] += 1
        dfRGI = pd.DataFrame(lRGI)
        dfRGI.to_excel(f"{WGS.seq_name}.xlsx")

    elif prgArgs.csvfile :
        logger.info(f"=[PROCESS] Multiple - {prgArgs.csvfile} [SEQ_VALID>0]")
        SEQ = pd.read_csv(prgArgs.csvfile)
        lstData = []
        for idx,row in SEQ.iterrows():
            
            if row['SEQ_VALID'] > 0 :
                logger.info(f"=[PROCESS] RGI - {row['ORGBATCH_ID']} {row['SEQRUN_ID']}")
                WGS = WGS_RDM(djDir['rdmDir_WGS'],row['ORGBATCH_ID'],row['SEQRUN_ID'])
                lstData += get_RGI_Info(WGS.fasta_dir,WGS.orgbatch_id, WGS.run_id)
                nProc['Processed'] += 1

        dfRGI = pd.DataFrame(lstData)
        dfRGI.to_excel(f"RGI_Data.xlsx")

    print(f"[WGS-Assembly] {nProc} ")


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading WGS Assembly to adjCOADD from RDM")

    prgParser = configargparse.ArgumentParser()
    prgParser.add_argument("-o","--orgbatch", default=None,required=False, dest="orgbatchid", action='store', help="OrgBatch_ID")
    prgParser.add_argument("-r","--runid", default=None,required=False, dest="runid", action='store', help="Seq RunID")

    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")

    prgParser.add_argument("-c","--csvfile", default=None,required=False, dest="csvfile", action='store', help="CSVFile")
    prgParser.add_argument("-e","--excel", default=None,required=False, dest="excelfile", action='store', help="ExcelFile")
    prgParser.add_argument("-s","--sheetname", default=None,required=False, dest="sheetname", action='store', help="ExcelSheet")

    prgParser.add_argument("-p","--process", default=None,required=False, dest="process", action='store', help="List of Processes [,] ",
                type=lambda s: [item for item in s.split(',')])

    #prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    #prgParser.add_argument("-l",default=None,required=True, dest="library", action='store', help="Library")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("--config",type=Path,is_config_file=True,help="Path to a configuration file ",)


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