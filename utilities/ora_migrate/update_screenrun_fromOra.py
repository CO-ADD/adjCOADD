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



def get_oraScreenRun():
    from oraCastDB.oraCastDB import openCastDB

    renameCol = {
        "run_type":   "run_project",
    }

    replaceValues = {
      'run_status':{'Completed Closed':'Completed', 
                    'Planned':'Planned', 
                    'Analysis Completed':'Reviewed', 
                    'Assay Uploaded':'Started', 
                    'Assay Completed':'Started', 
                    'Scheduled':'Planned', 
                    'Ongoing':'Ongoing'
                    },
      'screen_type':{'PS':'PSR', 
                     'HC':'HCR',
                     'HV':'HVR',
                     'SYN':'SYN',
                     'QC':'QC',
                     'Study':'STUDY'},
      'run_project':{'CooperGrp':'IMB',
                     'Internal':'IMB',
                     'CO-ADD':'CO-ADD',
                     'External':'Collab',
                     'Collaboration':'Collab',},
    }

    prjSQL = """
     Select run_id, run_name, run_conditions, run_issues, 
            run_type, screen_type, run_status, 
            project_lst, stock_format, cpoz_id, pub_status
     From Project
    """
    CastDB = openCastDB()
    logger.info(f"[ScreenRuns] ... ")
    runDF = pd.DataFrame(CastDB.get_dict_list(prjSQL))
    nTotal = len(runDF)
    logger.info(f"[ScreenRuns] {nTotal} ")
    CastDB.close()

    logger.info(f"DF - Rename Columns {len(renameCol)}")
    runDF.rename(columns=renameCol, inplace=True)

    logger.info(f"DF - Replace Values {len(replaceValues)}")
    for k in replaceValues:
        runDF[k].replace(replaceValues[k],inplace=True)
    return(runDF)


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_Fields_fromDict
    from dscreen.models import Screen_Run

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "RunID" :

        runDF = get_oraScreenRun()
        OutFile = f"UpdateScreenRun_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
        
        cpyFields = ['run_id','run_name',
                    'run_conditions','run_issues','run_project',
                    #'project_lst', 'stock_format', 'cpoz_id', 'pub_status',
                    ]
        arrayFields = {}
        dictFields = ['run_status','screen_type',]


        outNumbers = {'Proc':0,'New':0,'Upload':0}
        outDict = []    
        for idx,row in tqdm(runDF.iterrows(), total=runDF.shape[0]):
            new_entry = False
            outNumbers['Proc'] += 1
            
            # cvPrj = Convert_ProjectID.get(row['ora_project_id'])
            # if not cvPrj:
            #     cvPrj = Convert_ProjectID.new_COADD_Project_ID(row['ora_project_id'])

            # if cvPrj:
                #print(f"{row['ora_project_id']} {cvPrj.project_id}")
            djRun = Screen_Run.get(cvPrj.project_id)
            if djRun is None:
                djRun = Screen_Run()
                djRun.run_id = row['run_id']
                new_entry = True
                outNumbers['New'] += 1
            else:
                row['Issue'] = f"Exists"

            
            # set_dictFields(djPrj,row,cpyFields)
            # set_arrayFields(djPrj,row,arrayFields)
            # set_Dictionaries(djPrj,row,dictFields)

            validStatus = True
            validStatus = set_Fields_fromDict(djRun,row,
                                              FieldList=cpyFields, 
                                              ArrayDict=arrayFields, 
                                              DictList=dictFields)
            # djPrj.clean_Fields()
            # validDict = djPrj.validate()

            # if validDict:
            #     validStatus = False
            #     for k in validDict:
            #         print('Warning',k,validDict[k],'-')
            #     outDict.append(row)

            if validStatus:
                if prgArgs.upload:
                    if new_entry or prgArgs.overwrite:
                        outNumbers['Upload'] += 1
                        djRun.save()
                        
        print(f"[oraScreenRun] : {outNumbers}")
        if len(outDict) > 0:
            print(f"Writing {len(outDict)} Issues: {OutFile}")
            outDF = pd.DataFrame(outDict)
            outDF.to_excel(OutFile)
        else:
            print(f"No Issues")
    
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
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
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
