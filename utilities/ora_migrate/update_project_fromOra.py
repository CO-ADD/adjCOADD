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

def get_oraProject():
    from oraCastDB.oraCastDB import openCastDB

    renameCol = {
        "project_status":   "process_status",
        "group_id":         "ora_group_id",
        "organisation":     "ora_organisation",
        "report_ps_date":   "ora_psreport_date",
        "report_hc_date":   "ora_hcreport_date",
        "report_hv_date":   "ora_hvreport_date",
        "project_id":       "ora_project_id",
        "country":          "ora_country",
    }

    replaceValues = {
      'provided_container' : {'Tubes':'Tube'} 
    }

    prjSQL = """
     Select project_id, project_name, project_status, project_type, library_name,
            coadd_id, cpoz_id,
            received, completed,
            compound_comment, compound_status,
            group_id, contact_a_id, contact_b_id,
            organisation, country,
            report_comment, report_hc_date, report_hv_date, report_ps_date, report_status,
            screen_comment,
            screen_conc, screen_conc_unit, screen_status,
            stock_comment, stock_conc, stock_conc_unit, stock_container, stock_status,
            data_comment, data_status,
            antimicro_status,
            project_action, project_comment,
            provided_comment, provided_container,
            pub_date, pub_status
     From Project
    -- Where Organism_Name like 'Klebsiella%'
    """
    CastDB = openCastDB()
    logger.info(f"[Projects] ... ")
    prjDF = pd.DataFrame(CastDB.get_dict_list(prjSQL))
    nTotal = len(prjDF)
    logger.info(f"[Projects] {nTotal} ")
    CastDB.close()

    prjDF.rename(columns=renameCol, inplace=True)
    for k in replaceValues:
        prjDF[k].replace(to_replace=replaceValues[k],inplace=True)
    return(prjDF)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from dsample.models import Project
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------
    runTables = ["CompoundID","ProjectID"]

    if prgArgs.table == "ProjectID" :

        prjDF = get_oraProject()

        ignoreFields = ['COUNTRY','STATUS','PROJECT_ACTION','SCREEN_CONC','SCREEN_CONC_UNIT','COADD_ID','ANTIMICRO_STATUS']
        cpyFields = ['project_name','project_comment',
                    'provided_comment','stock_container','cpoz_id','process_status',
                    'received','completed',
                    'screen_comment','report_comment',
                    'compound_comment','data_comment',
                    'stock_comment', 'stock_conc',
                    'pub_name', 'pub_date',		
                    'ora_group_id','ora_organisation','ora_psreport_date','ora_hcreport_date','ora_hvreport_date','ora_project_id'
                    ]
        arrayFields = {'screen_status': 'screen_status',
                        'report_status': 'report_status',
                        'compound_status': 'compound_status',
                        'data_status': 'data_status',
                        'stock_status': 'stock_status',
                        'pub_status': 'pub_status',
                        'ora_contact_ids':['CONTACT_A_ID','CONTACT_B_ID']
                    }



        # ExcelFile = prgArgs.file
        # SheetName = "ProjectID"
        # print(f"[Reading Excel] {ExcelFile} {SheetName} ")
        # prjDF = pd.read_excel(ExcelFile,sheet_name=SheetName)
        # prjDF.project_name = prjDF.project_name.fillna('')
        
        # for idx,row in tqdm(prjDF.iterrows(), total=prjDF.shape[0]):
        #     djE = Convert_ProjectID.get(row['ora_project_id'])
        #     if djE is None:
        #         djE = Convert_ProjectID()
        #         djE.ora_project_id = row['ora_project_id']
        #         djE.project_id = row['project_id']
        #         djE.project_name = row['project_name']
        #         if prgArgs.upload:
        #             djE.clean_Fields()
        #             djE.save()
            # else:
            #     print(f"[Exists already] {row['ora_project_id']} {row['project_id']} ")

    
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
        djDir = "D:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
