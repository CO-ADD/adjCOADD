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
logName = "Sum_CmpBatch_Inhib"
logFileName = os.path.join("log",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------

def get_strList(x): 
    return ";".join(x) 

def get_nAct(x):
    return len([a for a in x if a == 'A']) 

def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dplate.models import Labware, TestPlate, TestWell
    from dsample.models import COADD_Compound, Compound_Batch
    from dsummary.utils.upd_sum_cmpbatch import sum_structure_sc
    from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
    from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp
    from adjcoadd.constants import COMPOUND_SEP
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # AssayData MIC -------------------------------------------------------------
    if prgArgs.table == 'Sum_Structure_SC':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"{prgArgs.table}_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        qrySources = ['COADD']

        if int(prgArgs.test) > 0:
            twStr = TestWell.objects.filter(n_cmpbatches = 1, mpbatch_id__structure_id__isnull = False, 
                                            plate_id__result_type = 'Inhibition').values_list('cmpbatch_id__structure_id').distinct()[:int(prgArgs.test)]
        else:
            twStr = TestWell.objects.filter(n_cmpbatches = 1, mpbatch_id__structure_id__isnull = False, 
                                            plate_id__result_type = 'Inhibition').values_list('cmpbatch_id__structure_id').distinct()

        logger.info(f" [Sum Structure SC] TestWell: {twStr.count()}  ")

        for sid in tqdm(twStr, desc='[Structures]'):
            #print(cmps)
            _numbers,_outdict  = sum_structure_sc(sid[0],upload=prgArgs.upload,overwrite=prgArgs.overwrite,appuser=prgArgs.appuser)

            if _outdict:
                OutDict = OutDict + _outdict
            for k in OutNumbers.keys():
                OutNumbers[k] += _numbers[k]

        #         if _outdict:
        #             OutDict = OutDict + _outdict
        #         for k in OutNumbers.keys():
        #             OutNumbers[k] += _numbers[k]

        #     if len(OutDict) > 0:
        #         logger.info(f"Writing Issues: {OutFile}")
        #         outDF = pd.DataFrame(OutDict)
        #         outDF.to_excel(OutFile)
        #     else:
        #         logger.info(f"No Issues")


        if len(OutDict) > 0:
            logger.info(f"Writing Issues: {OutFile}")
            outDF = pd.DataFrame(OutDict)
            outDF.to_excel(OutFile)
        else:
            logger.info(f"No Issues")

        logger.info(f"{OutName} {OutNumbers}")


        # for cmp in tqdm(qryCmp, desc='[Compounds]'):
        #     validStatus = True
        #     CmpBatchID = cmp[0]
        #     qryCmpBatchID = [CmpBatchID]
        #     qryNCmpBatches = len(qryCmpBatchID)

        #     qryTW = TestWell.objects.filter(cmpbatch_lst__contains = [CmpBatchID], 
        #                                     n_cmpbatches = qryNCmpBatches, 
        #                                     plate_id__result_type = 'Inhibition',
        #                                     is_valid = True,
        #                                     plate_id__plate_quality = 'Valid'
        #                                    ).values(
        #                                        'plate_id','well_id','plate_id__result_type','plate_id__assay_id',
        #                                        'inhibition','mscore','act_type'
        #                                          )
        #     if qryTW.exists():
        #         dfSC = pd.DataFrame(qryTW).assign(cmpbatch_id=CmpBatchID)
        #         dfSC.columns = ['plate_id','well_id','result_type','assay_id','inhibition','mscore','act_type','cmpbatch_id']

                #print(dfSC)

                # pivDF = dfSC.groupby(['assay_id']).agg({'inhibition': ['mean','max','min','std'],
                #                                     'mscore': ['mean','size'],
                #                                     'act_type': [lambda x: ";".join(x), lambda x: len(x) if (x == "A").any() else 0 ],
                #                                      })
                # pivDF = dfSC.groupby(['assay_id']).agg({'inhibition': ['mean','max','min','std'],
                #                                     'mscore': ['mean','size'],
                #                                     'act_type': [get_strList, get_nAct ],
                #                                     })
                
                # #print( pivDF.columns)


                # for idx,row in pivDF.iterrows():
                #     OutNumbers['Processed'] += 1
                #     #print(idx," --> ", row.to_dict())
                #     NewEntry = False
                #     djSum = Summary_CmpBatch_Inhib.get(qryCmpBatchID,idx,Exact=True,verbose=0)
                #     if djSum is None:
                #         djSum = Summary_CmpBatch_Inhib()
                #         djSum.cmpbatch_lst = qryCmpBatchID
                #         djSum.n_cmpbatches = len(qryCmpBatchID)
                #         djSum.assay_id = idx
                #         NewEntry = True
                #     djSum.act_types = row[ ('act_type','get_strList')]

                #     djSum.n_assays = row[('mscore','size')]
                #     djSum.n_actives = row[ ('act_type','get_nAct')]
                #     #djSum.act_score_ave =

                #     djSum.inhibition_ave = row[('inhibition','mean')]
                #     djSum.inhibition_std = row[('inhibition','std')]
                #     djSum.inhibition_min = row[('inhibition','min')]
                #     djSum.inhibition_max = row[('inhibition','max')]
                #     djSum.mscore_ave = row[('mscore','mean')]


                #     djSum.clean_Fields()
                #     validDict = djSum.validate()
                #     if validDict:
                #         validStatus = False
                #         # for k in validDict:
                #         #     print('Warning',k,validDict[k],'-')
                #         row.update(validDict)

                #     if validStatus:
                #         if prgArgs.upload:
                #             if NewEntry or prgArgs.overwrite:
                #                 #djSum.chk_migration = 0
                #                 OutNumbers['Upload Entries'] += 1
                #                 djSum.save(user=prgArgs.appuser)
                                

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