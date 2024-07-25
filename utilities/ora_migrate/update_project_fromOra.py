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
        "library_name":     "pub_name",
    }

    replaceValues = {
      'provided_container':{'Tubes':'Tube', 'Vials':'Vial', 'Plates':'Plate'},
      'project_type':{'Internal':'Reference', 'Project':'GrpProject', 'Agreement':'Contract'},
      'stock_conc_unit':{'mg':'mg/mL'},
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

    logger.info(f"DF - Rename Columns {len(renameCol)}")
    prjDF.rename(columns=renameCol, inplace=True)

    logger.info(f"DF - Replace Values {len(replaceValues)}")
    for k in replaceValues:
        prjDF[k].replace(replaceValues[k],inplace=True)
    return(prjDF)


#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from dsample.models import Project
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "ProjectID" :

        prjDF = get_oraProject()
        OutFile = f"UpdateProject_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"
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
                        'ora_contact_ids':['contact_a_id','contact_b_id']
                    }
        dictFields = ['project_type','provided_container','stock_conc_unit',]


        outNumbers = {'Proc':0,'New':0,'Upload':0}
        outDict = []    
        for idx,row in tqdm(prjDF.iterrows(), total=prjDF.shape[0]):
            new_entry = False
            outNumbers['Proc'] += 1
            cvPrj = Convert_ProjectID.get(row['ora_project_id'])

            if not cvPrj:
                cvPrj = Convert_ProjectID.new_COADD_Project_ID(row['ora_project_id'])

            if cvPrj:
                #print(f"{row['ora_project_id']} {cvPrj.project_id}")
                djPrj = Project.get(cvPrj.project_id)
                if djPrj is None:
                    djPrj = Project()
                    djPrj.project_id = cvPrj.project_id
                    new_entry = True
                    outNumbers['New'] += 1
                else:
                    row['Issue'] = f"Exists"

                set_dictFields(djPrj,row,cpyFields)
                set_arrayFields(djPrj,row,arrayFields)
                set_Dictionaries(djPrj,row,dictFields)

                validStatus = True
                djPrj.clean_Fields()
                validDict = djPrj.validate()

                if validDict:
                    validStatus = False
                    for k in validDict:
                        print('Warning',k,validDict[k],'-')
                    outDict.append(row)

                if validStatus:
                    if prgArgs.upload:
                        if new_entry or prgArgs.overwrite:
                            outNumbers['Upload'] += 1
                            djPrj.save()
            else:
                row['Issue'] = f"ConvProject not found"
                print(f"[oraProject] ConvProject {row['ora_project_id']} not found")
                outDict.append(row)
        print(f"[oraProject] : {outNumbers}")
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
