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
logName = "Sum_Structure_DoseResp"
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

    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from dsummary.utils.upd_sum_cmpbatch import sum_structure_dr
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from adjcoadd.constants import COMPOUND_SEP

   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'Sum_Structure_DR':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"{prgArgs.table}_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        # Get all Distinct CmpBatch_Lst
        if prgArgs.structureid:
            strDict = {prgArgs.structureid:prgArgs.structureid}
        else:
            if int(prgArgs.test) > 0:
                micStruct = AssayData_MIC.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()[:int(prgArgs.test)]
                cc50Struct = AssayData_CC50.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()[:int(prgArgs.test)]
                hc50Struct = AssayData_HC50.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()[:int(prgArgs.test)]
            else:
                micStruct = AssayData_MIC.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()
                cc50Struct = AssayData_CC50.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()
                hc50Struct = AssayData_HC50.objects.filter(n_cmpbatches = 1, 
                                                        cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()

            logger.info(f" [Sum Structure DR] MIC: {micStruct.count()} + CC50: {cc50Struct.count()} + HC50: {hc50Struct.count()} ")

            # Get unique ID
            strDict = {}
            for s in micStruct:
                if s[0] not in strDict:
                    strDict[s[0]] = s[0]
            for s in cc50Struct:
                if s[0] not in strDict:
                    strDict[s[0]] = s[0]
            for s in hc50Struct:
                if s[0] not in strDict:
                    strDict[s[0]] = s[0]
            logger.info(f" [Sum Structure DR] Structures: {len(strDict)} ")

        for sid in tqdm(strDict.keys(), desc='[CmpBatcheLsts]'):
            #print(cmpDict[cmps])
            _numbers,_outdict  = sum_structure_dr(sid,upload=prgArgs.upload,overwrite=prgArgs.overwrite,appuser=prgArgs.appuser)

            if _outdict:
                OutDict = OutDict + _outdict
            for k in OutNumbers.keys():
                OutNumbers[k] += _numbers[k]

        if len(OutDict) > 0:
            logger.info(f"Writing Issues: {OutFile}")
            outDF = pd.DataFrame(OutDict)
            outDF.to_excel(OutFile)
        else:
            logger.info(f"No Issues")

        logger.info(f"{OutName} {OutNumbers}")


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
    prgParser.add_argument("-s",default=None,required=False, dest="structureid", action='store', help="Single File to parse")

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