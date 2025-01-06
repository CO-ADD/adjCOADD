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
logName = "Sum_CmpBatch_DoseResp"
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

    from django.db.models import Q
    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries, set_fkeyFields, set_arrayDictionaries
    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from dsummary.utils.upd_sum_cmpbatch import sum_cmpbatch_dr
    from ddrug.utils.bio_data import DR_Range
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from dorganism.models import Organism_Batch
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
#    from dorganism.utils.utils  import reformat_OrganismID, reformat_OrgBatchID
#    from update_utils import convert_castdb_compoundid_from_ora
    from adjcoadd.constants import COMPOUND_SEP

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Agg Funvtions  -------------------------------------------------------------
    def get_strList(x, maxN = 10):
        if len(x) > maxN:
            _v, _c = np.unique(x, return_counts=True)
            _a = []
            for _i in range(len(_v)):
                _a.append(f"{_v[_i]} ({_c[_i]})")
            return "; ".join(_a)     
        return "; ".join(x) 

    def get_strList_unique(x):
        return ";".join(set(x)) 

    def get_nAct(x):
        return len([a for a in x if a == 'A']) 

    def get_DR_Range(x):
        return DR_Range(x)
    
   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'Sum_CmpBatch_DR':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"{prgArgs.table}_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        # Get all Distinct CmpBatch_Lst
        if int(prgArgs.test) > 0:
            micCmp = AssayData_MIC.objects.all().values_list('cmpbatch_lst').distinct()[:int(prgArgs.test)]
            cc50Cmp = AssayData_CC50.objects.all().values_list('cmpbatch_lst').distinct()[:int(prgArgs.test)]
            hc50Cmp = AssayData_HC50.objects.all().values_list('cmpbatch_lst').distinct()[:int(prgArgs.test)]
        else:
            micCmp = AssayData_MIC.objects.all().values_list('cmpbatch_lst').distinct()
            cc50Cmp = AssayData_CC50.objects.all().values_list('cmpbatch_lst').distinct()
            hc50Cmp = AssayData_HC50.objects.all().values_list('cmpbatch_lst').distinct()


        logger.info(f" [Sum CmpBatch DR] MIC: {micCmp.count()} + CC50: {cc50Cmp.count()} + HC50: {hc50Cmp.count()} ")

        # Distinct CmpBatch_Lst
        cmpDict = {}
        for c in micCmp:
            cc = COMPOUND_SEP.join([str(x) for x in c[0] if x != ""])
            if cc not in cmpDict:
                cmpDict[cc] = c[0]
        for c in cc50Cmp:
            cc = COMPOUND_SEP.join([str(x) for x in c[0] if x != ""])
            if cc not in cmpDict:
                cmpDict[cc] = c[0]
        for s in hc50Cmp:
            cc = COMPOUND_SEP.join([str(x) for x in c[0] if x != ""])
            if cc not in cmpDict:
                cmpDict[cc] = c[0]

        logger.info(f" [Sum CmpBatch DR] CmpBatcheLsts: {len(cmpDict)} ")

        for cmps in tqdm(cmpDict.keys(), desc='[CmpBatcheLsts]'):
            #print(cmpDict[cmps])
            _numbers,_outdict  = sum_cmpbatch_dr(cmpDict[cmps],upload=prgArgs.upload,overwrite=prgArgs.overwrite,appuser=prgArgs.appuser)

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