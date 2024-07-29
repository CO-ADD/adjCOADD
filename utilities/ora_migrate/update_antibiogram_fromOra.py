#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from oraCastDB.oraCastDB import openCastDB

from tqdm import tqdm
# from zUtils import zData

import django

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_ConvertID"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)

#-----------------------------------------------------------------------------


def get_MIC_COADD(test=0):
    
    renameCol = {
        "test_dye":      "dye",
        "test_additive":      "additive",
        "media_id":       "media",
        "material":    "plate_material",
    }

    replaceValues = {
      'reg_conc_unit':{'ug/ul':'mg/mL','mg/ml':'mg/mL'},
    }

    lstRunID =['AntiBio_R001','AntiBio_R002','AntiBio_R003',
               'AntiBio_R004','AntiBio_R005','AntiBio_R006',
               'AntiBio_R007','AntiBio_R008','AntiBio_R009',
               'HCR00132', 
               'HVR00049','HVR00051','HVR00052','HVR00057',
               'PMC043_R02','PMC043_R03',
               'PMC045_FDB1',
               'PMC045_R01','PMC045_R02','PMC045_R03',
               'PMC045_R04','PMC045_R05','PMC045_R06']

    micSQL = """
            Select mic.TestPlate_ID, mic.TestWell_ID, 
                mic.Compound1_ID, mic.Compound2_ID, mic.Compound3_ID, mic.Compound4_ID,
                c1.Compound_Name Compound1_Name, c2.Compound_Name Compound2_Name, 
                c3.Compound_Name Compound3_Name, c4.Compound_Name Compound4_Name,
                mic.AssayType_Id, 
                tp.Test_Strain, tp.Plate_Size, tp.Media_id, lw.Material, tp.Test_Dye, tp.Test_Additive, tp.Test_Date,
                mic.MIC, mic.MIC_Unit, mic.DMax, mic.Data_Quality, mic.Run_ID 
            From AssayData_MIC mic
            Left Join Compound c1 on c1.Compound_ID = mic.Compound1_ID
            Left Join Compound c2 on c2.Compound_ID = mic.Compound2_ID
            Left Join Compound c3 on c3.Compound_ID = mic.Compound3_ID
            Left Join Compound c4 on c4.Compound_ID = mic.Compound4_ID
            Left Join Testplate tp on tp.Plate_id = mic.Testplate_ID
            Left Join Labware lw on lw.Labware_id = tp.Labware_id
            Where c1.project_id in ('PC001','PC002') and mic.n_compounds < 3
                """
    
    _sqlRun = ""            
    for r in lstRunID:
        _sqlRun += f"'{r}',"
    micSQL += f" and mic.Run_ID in ({_sqlRun[:-1]})"


    if test>0:
        micSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[MIC-COADD] ... ")
    micDF = pd.DataFrame(CastDB.get_dict_list(micSQL))
    #micLst = CastDB.get_dict_list(micSQL)
    nTotal = len(micDF)
    logger.info(f"[MIC-COADD] {nTotal} ")
    CastDB.close()


    logger.info(f"DF - Rename Columns {len(renameCol)}")
    micDF.rename(columns=renameCol, inplace=True)

    # logger.info(f"DF - Replace Values {len(replaceValues)}")
    # for k in replaceValues:
    #     cmpDF[k].replace(replaceValues[k],inplace=True)

    return(micDF)

# #-----------------------------------------------------------------------------------
# def update_MICCOADD_ora(RunID,upload=False,uploaduser=None,OutputN=100):
# #-----------------------------------------------------------------------------------
#     if nTotal>0:
#         vLog = validation_log.Validation_Log('MIC-Pub')
#         nTime  = zData.Timer(nTotal)
#         nProcessed = 0

#         # check user
#         appuser = None
#         if uploaduser:
#             appuser = ApplicationUser.get(uploaduser)

#         nProc = {}
#         nProc['Saved'] = 0
#         nProc['notValid'] = 0

#         for mic in micLst:
#             mic['ORGBATCH_ID'] = reformat_OrgBatchID(mic['TEST_STRAIN'])
#             if mic['COMPOUND2_NAME']:
#                 mic['DRUG_NAME'] = mic['COMPOUND_NAME'] +'|'+ mic['COMPOUND2_NAME']
#                 if mic['DR2_VALUE']:
#                     mic['MIC'] = mic['DR'] +'|'+ str(mic['DR2_VALUE'])
#                 else:
#                     mic['MIC'] = mic['DR'] 
#                 if mic['DR2_UNIT']:
#                     mic['MIC_UNIT'] = mic['DR_UNIT'] +'|'+ mic['DR2_UNIT']
#                 else:
#                     mic['MIC_UNIT'] = mic['DR_UNIT']
#             else:
#                 mic['DRUG_NAME'] = mic['COMPOUND_NAME']
#                 mic['MIC'] = mic['DR']
#                 mic['MIC_UNIT'] = mic['DR_UNIT']
            
#             mic['MEDIA'] = None
#             mic['PLATE_SIZE'] = mic['PLATE_SIZE'].replace('w','')
#             mic['PLATE_MATERIAL'] = mic['MATERIAL']
#             mic['DYE'] = mic['TEST_DYE']
#             mic['ADDITIVE'] = mic['TEST_ADDITIVE']

#             djMIC = imp_MICCOADD_fromDict(mic,vLog)
#             djMIC.clean_Fields()
#             validDict = djMIC.validate()

#             if validDict:
#                 logger.info(f"{mic['ORGBATCH_ID']} {mic['DRUG_NAME']} {mic['RUN_ID']} {validDict} ")

#             # --- Upload ---------------------------------------------------------
#             nProcessed = nProcessed + 1
#             if djMIC.VALID_STATUS:
#                 if upload:
#                     if nProcessed%OutputN == 0:
#                         eTime,sTime = nTime.remains(nProcessed)
#                         logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} -> {djMIC} ")
#                     djMIC.save(user=appuser)
#                     nProc['Saved'] = nProc['Saved'] + 1
#                 else:
#                     if nProcessed%OutputN == 0:
#                         eTime,sTime = nTime.remains(nProcessed)
#                         logger.info(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >r {djMIC} ")
#             else:
#                 nProc['notValid'] = nProc['notValid'] + 1
#         eTime,sTime = nTime.remains(nProcessed)
#         logger.info(f"[MIC-COADD] [{nTotal-(nProc['Saved']+nProc['notValid'])}] {nTotal} {nProc}")
#     else:
#         logger.info(f"[MIC-COADD] [0 No records found]")


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "regChem_01eRegChem_COADD"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
        #handlers=[logging.StreamHandler()],
        level=logLevel)
    #-----------------------------------------------------------------------------


    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.data import join_lst
    from apputil.utils.set_data import set_Fields_fromDict
    from dorganism.utils.utils import reformat_OrgBatchID

    from ddrug.models import Drug, MIC_COADD
    from dorganism.models import Organism_Batch
    from dscreen.models import Screen_Run
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "Antibiogram" :

        micDF = get_MIC_COADD(int(prgArgs.test))
        print("--> oraMIC ---------------------------------------------------------")

        print(micDF.columns)
        print("-------------------------------------------------------------------------")
        OutFile = f"UpdateMIC_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"


        outNumbers = {'Proc':0,'Upload':0, 'New Data':0,'Invalid Data':0,'Missing Data':0}
        outDict = []    
        for idx,row in tqdm(micDF.iterrows(), total=micDF.shape[0]):
            outNumbers['Proc'] += 1
            if 'Valid' in row['data_quality'] or 'Retest' in row['data_quality']:
                new_data = False
                

                drugName = join_lst([row['compound1_name'],row['compound2_name'],row['compound3_name'],row['compound4_name']],sep='|')
                orgbatchid = reformat_OrgBatchID(row['test_strain'])
                runid = row['run_id']
                
                validStatus = True
                djDrug = Drug.get(drugName)
                if djDrug is None:
                    validStatus = False
                    print(f" [Drug] not found '{drugName}'")

                djOrgBatch = Organism_Batch.get(orgbatchid) 
                if djOrgBatch is None:
                    validStatus = False
                    print(f" [OrgBatchID] not found '{orgbatchid}'")

                djRun = Screen_Run.get(runid)
                if djRun is None:
                    validStatus = False
                    print(f" [RunID] not found '{runid}'")

                if validStatus:
                    # 
                    # get(cls,OrgBatchID,DrugID,TestPlateID,TestWellID,verbose=0)
                    
                    djMIC = MIC_COADD.get(orgbatchid,str(djDrug),row['testplate_id'],row['testwell_id'])
                    if djMIC is None:
                        djMIC = MIC_COADD()
                        djMIC.testplate_id = row['testplate_id']
                        djMIC.testwell_id = row['testwell_id']
                        new_data = True
                        outNumbers['New Data'] += 1
                        #print(f" [MIC] new '{row['testplate_id']}:{row['testwell_id']}'")
                    djMIC.run_id = djRun
                    djMIC.orgbatch_id = djOrgBatch
                    djMIC.drug_id = djDrug
                    row['plate_size'] = row['plate_size'].replace('w','')
                    row['mic_type'] = 'BMD'
                    
                    validStatus = set_Fields_fromDict(djMIC,row,
                                                    FieldList=['mic','mic_unit','media','dye','additive','test_date'], 
                                                    ArrayDict={}, 
                                                    DictList=['plate_size','plate_material','mic_type'])

                    if validStatus:
                        if prgArgs.upload:
                            if new_data or prgArgs.overwrite:
                                outNumbers['Upload'] += 1
                                djMIC.save()
                else:
                    outNumbers['Missing Data'] += 1
                    print(f" [Antibiogram] invalid {drugName} {orgbatchid} {runid}")
            else:
                outNumbers['Invalid Data'] += 1
                 
        print(f"[Antibiogram] : {outNumbers}")



#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of rows to test")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--config",default='Local',required=False, dest="config", action='store', help="Configuration [Meran/Laptop/Work]")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
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
