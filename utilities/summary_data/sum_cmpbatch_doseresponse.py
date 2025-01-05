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
    from dsummary.utils.upd_sum_cmpbatch import sum_cmpbatch_doseresponse
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
    if prgArgs.table == 'Summary_CmpBatch_DoseResponse':

        OutName = f"[{prgArgs.table}]"
        OutDict = []
        OutFile = f"{prgArgs.table}_{logTime:%Y%m%d_%H%M%S}.xlsx"
        OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}

        qrySources = ['COADD']

        if int(prgArgs.test) > 0:
            qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')[:int(prgArgs.test)]
            qryCmpBatch = AssayData_MIC.objects.order_by().values_list('cmpbatch_lst').distinct()[:int(prgArgs.test)]
        else:
            #qryCmp = Compound_Batch.objects.filter(batch_source__in=qrySources).values_list('cmpbatch_id')
            qryCmpBatch = AssayData_MIC.objects.order_by().values_list('cmpbatch_lst').distinct()
        
        for cmp in tqdm(qryCmpBatch, desc='[Compounds]'):
            validStatus = True
            qryCmpBatchLst = cmp[0]

            _numbers,_outdict  = sum_cmpbatch_doseresponse(cmp[0],upload=prgArgs.upload,overwrite=prgArgs.overwrite,appuser=prgArgs.appuser)

            OutDict = OutDict + _outdict
            for k in OutNumbers.keys():
                OutNumbers[k] += _numbers[k]

            # CmpBatchs = COMPOUND_SEP.join(qryCmpBatchLst)
            # qryNCmpBatches = len(qryCmpBatchLst)

            # # Sum_Cmpd ---------------------------------------------------------------
            # djCmpd = Summary_CmpBatch.get(qryCmpBatchLst,verbose=0)
            # if djCmpd is None:
            #     djCmpd = Summary_CmpBatch()
            #     djCmpd.cmpbatch_lst = qryCmpBatchLst
            #     djCmpd.n_cmpbatches = len(qryCmpBatchLst)

            # djCmpd.dr_n_assayids = 0
            # djCmpd.dr_n_actives = 0

            # # AssayData  ---------------------------------------------------------------
            # qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
            #                                 cmpbatch_lst__contains = qryCmpBatchLst, 
            #                                 n_cmpbatches = qryNCmpBatches, 
            #                                 testplate_id__result_type = 'MIC',
            #                                 testplate_id__plate_quality = 'Valid'                                            
            #                                ).values(
            #                                    'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
            #                                    'mic','mic_unit','act_type','act_score','pscore',
            #                                    'inhibit_max'
            #                                      )
            # if qryMIC.exists():
            #     dfDR = pd.DataFrame(qryMIC).assign(cmpbatchs=CmpBatchs)
            #     dfDR.columns = ['plate_id','well_id','result_type','assay_id',
            #                     'mic','mic_unit','act_type','act_score','pscore',
            #                     'inhibit_max',
            #                     'cmpbatchs']

            #     #print(dfDR)
            #     pivDF = dfDR.groupby(['assay_id']).agg({'mic': [DR_Range],
            #                                         'inhibit_max': ['mean'],
            #                                         'pscore': ['mean'],
            #                                         'act_score': ['mean'],
            #                                         'act_type': [get_strList, get_nAct ],
            #                                         'mic_unit': [get_strList, get_strList_unique ],
            #                                         })
            #     #print( pivDF.columns)
            #     for idx,row in pivDF.iterrows():
            #         OutNumbers['Processed'] += 1
            #         djCmpd.dr_n_assayids += 1

            #         NewEntry = False
            #         djSum = Summary_CmpBatch_Doseresp.get(qryCmpBatchLst,idx,Exact=True,verbose=0)
            #         if djSum is None:
            #             djSum = Summary_CmpBatch_Doseresp()
            #             djSum.cmpbatch_lst = qryCmpBatchLst
            #             djSum.n_cmpbatches = len(qryCmpBatchLst)
            #             djSum.assay_id = idx
            #             NewEntry = True
            #         djSum.act_types = row[('act_type','get_strList')]
            #         djSum.n_actives = row[('act_type','get_nAct')]
            #         djSum.act_score_ave = row[('act_score','mean')]
            #         djSum.inhibit_max_ave = row[('inhibit_max','mean')]
            #         djSum.pscore = row[('pscore','mean')]
            #         djSum.drval_type = 'MIC'

            #         if djSum.n_actives > 0:
            #             djCmpd.dr_n_actives += 1
            #         #print(f" {row[('mic','DR_Range')]}")

            #         djSum.drval_max    = row[('mic','DR_Range')]['Max']
            #         djSum.drval_min    = row[('mic','DR_Range')]['Min']
            #         djSum.drval_median = row[('mic','DR_Range')]['Median']
            #         djSum.n_assays = row[('mic','DR_Range')]['nDR']
            #         djSum.drval_unit   = row[('mic_unit','get_strList_unique')]

            #         djSum.clean_Fields()
            #         validDict = djSum.validate()
            #         if validDict:
            #             validStatus = False
            #             # for k in validDict:
            #             #     print('Warning',k,validDict[k],'-')
            #             row.update(validDict)
            #             OutDict.append(row)

            #         if validStatus:
            #             if prgArgs.upload:
            #                 if NewEntry or prgArgs.overwrite:
            #                     #djSum.chk_migration = 0
            #                     OutNumbers['Upload Entries'] += 1
            #                     djSum.save(user=prgArgs.appuser)

            # # Sum_Cmpd ---------------------------------------------------------------
            # # djCmpd.dr_assayid_lst  = 
            # # djCmpd.dr_actives_lst  = 
            # if prgArgs.upload:
            #     djCmpd.save()

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