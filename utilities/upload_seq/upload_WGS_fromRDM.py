#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

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

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from apputil.utils.data import listFolders
    from apputil.utils import validation_log
    
    from dgene.models import Gene,ID_Pub,ID_Sequence,WGS_FastQC,WGS_CheckM
    from dgene.utils.upload_gene import (get_RDM, split_BatchID_RunID, get_subdir,
                                        upload_Trim, upload_CheckM, upload_FastA, upload_AMR)
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

    if prgArgs.directory:
        RDM = get_RDM(prgArgs.directory)
        RDM['base'] = prgArgs.directory
        
        AssemblyBase = os.path.join(RDM['base'],RDM['assembly'])
        FastABase = os.path.join(RDM['base'],RDM['fasta'])

        vLog = validation_log.Validation_Log('WGS-Assembly')
        
        # for subDir in listFolders(AssemblyBase):
        #     zAssemblyFolder = os.path.join(AssemblyBase,subDir)
        #     for BatchRunID in listFolders(zAssemblyFolder):
        #         dirAss = os.path.join(zAssemblyFolder,f"{BatchRunID}")
        #         if os.path.exists(dirAss):

        #             #nProcessed = nProcessed + 1
        #             OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)
        #             if prgArgs.runid:
        #                 fProcess = prgArgs.runid == RunID
        #             else:
        #                 fProcess = True
                    
        #             if fProcess:
        #                 print(f"[WGS-Assembly] {OrgBatchID} {RunID} {dirAss} ")                       
        #                 upload_CheckM(OrgBatchID, RunID, dirAss, vLog, upload=prgArgs.upload,uploaduser=prgArgs.appuser)
                        

        Methods= ['AMR Finder']
        for subDir in listFolders(FastABase):
            zFastAFolder = os.path.join(FastABase,subDir)
            for BatchRunID in listFolders(zFastAFolder):
                dirFA = os.path.join(zFastAFolder,f"{BatchRunID}")
                if os.path.exists(dirFA):

                    OrgBatchID, RunID = split_BatchID_RunID(BatchRunID)

                    if prgArgs.runid:
                        fProcess = prgArgs.runid == RunID
                    else:
                        fProcess = True
                    
                    if fProcess:
                        print(f"[WGS-FastA] {OrgBatchID} {RunID} ")                       
                        upload_FastA(OrgBatchID, RunID, dirFA, vLog, upload=prgArgs.upload,uploaduser=prgArgs.appuser) 
                        upload_AMR(OrgBatchID, RunID, dirFA, vLog, Methods, upload=prgArgs.upload,uploaduser=prgArgs.appuser)
                                           
    print(f"[WGS-Assembly] {nProc} ")

    # if prgArgs.orgbatch and prgArgs.runid:
    #     dGene.update_WGSCOADD_Assembly_single(prgArgs.orgbatch,prgArgs.runid,upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    # else:   
    #     logger.info(f"[Upd_djCOADD] {prgArgs.table} from 02_Assembly {prgArgs.runid} [Upload: {prgArgs.upload}]")
    #     #dGene.update_WGSCOADD_Trim(upload=prgArgs.upload,uploaduser=prgArgs.appuser)
    #     dGene.update_WGSCOADD_Assembly(upload=prgArgs.upload,uploaduser=prgArgs.appuser)




#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading WGS Assembly to adjCOADD from RDM")
    #prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    #prgParser.add_argument("-l",default=None,required=True, dest="library", action='store', help="Library")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.config == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.config == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.config == 'Laptop':
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
