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

    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from adjcoadd.constants import COMPOUND_SEP

   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'Summary_Structure_DoseResponse':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"{prgArgs.table}_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}


        micStruct = AssayData_MIC.objects.filter(n_cmpbatches = 1, cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()
        cc50Struct = AssayData_CC50.objects.filter(n_cmpbatches = 1, cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()
        hc50Struct = AssayData_HC50.objects.filter(n_cmpbatches = 1, cmpbatch_id__structure_id__isnull = False ).values_list('cmpbatch_id__structure_id').distinct()

        print(f" MIC: {micStruct.count()} + CC50: {cc50Struct.count()} + HC50: {hc50Struct.count()} ")

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
        print(f" Structures: {len(strDict)} ")

        micCmp = AssayData_MIC.objects.all().values_list('cmpbatch_lst').distinct()
        cc50Cmp = AssayData_CC50.objects.all().values_list('cmpbatch_lst').distinct()
        hc50Cmp = AssayData_HC50.objects.all().values_list('cmpbatch_lst').distinct()

        print(f" MIC: {micCmp.count()} + CC50: {cc50Cmp.count()} + HC50: {hc50Cmp.count()} ")

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
        print(f" CmpBatches: {len(cmpDict)} ")



        #print(f" MIC: {len(list(micStruct))} + CC50: {len(list(cc50Struct))} + HC50: {len(list(hc50Struct))}")


        # if int(prgArgs.test) > 0:
        #     qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')[:int(prgArgs.test)]
        #     qryCmpBatch = AssayData_MIC.objects.order_by().values_list('cmpbatch_lst').distinct()[:int(prgArgs.test)]
        # else:
        #     #qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')
        #     qryCmpBatch = AssayData_MIC.objects.order_by().values_list('cmpbatch_lst').distinct()

        # for sid in tqdm(qryStruct, desc='[Structures]'):

        #     qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
        #                                     cmpbatch_id__structure_id = sid, 
        #                                     n_cmpbatches = 1, 
        #                                     testplate_id__result_type = 'MIC',
        #                                     testplate_id__plate_quality = 'Valid'                                            
        #                                     ).values(
        #                                         'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
        #                                         'mic','mic_unit','act_type','act_score','pscore',
        #                                         'inhibit_max'
        #                                             )


        # for row in qryStruct:
        #     print(row.cmpbatch_lst)

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