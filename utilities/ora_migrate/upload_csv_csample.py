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

import logging
#-----------------------------------------------------------------------------

def set_arrayFields_fromDict(djModel,rowDict, arrDict):
    for f in arrDict:
        #print("arrFields",f)
        if isinstance(arrDict[f],str):
            if pd.notnull(rowDict[arrDict[f]]):
                setattr(djModel,f,rowDict[arrDict[f]].split(";"))
                
        elif isinstance(arrDict[f],list):
            _list = []
            for l in arrDict[f]:
                if pd.notnull(rowDict[l]):
                    _list.append(rowDict[l])
            setattr(djModel,f,_list)
        
def set_fromDict(djModel,rowDict,setList):
    for e in setList:
        #print("fromDict",e)
        if pd.notnull(rowDict[e]):
            setattr(djModel,e,rowDict[e])

        
#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "D:/Code/zdjCode/adjCOADD"
    # uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    # if prgArgs.database == 'Work':
    #     djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    # elif prgArgs.database == 'WorkLinux':
    #     djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"

    # xlFiles = {
    #     'Application': "ApplicationData_v05.xlsx",
    #     'Drug': "DrugData_v04.xlsx",
    #     'MIC': "LMIC_Data_v06.xlsx",
    #     'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    # }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from dsample.models import Project, Convert_ProjectID, Sample, COADD_Compound, Convert_CompoundID
    from dchem.utils.mol_std import get_atomclass_list,list_metalatoms
    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "UploadOrgDB"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
    logLevel = logging.INFO 

    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="[%(name)-20s] %(message)s ",
        handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
        level=logLevel)
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    def set_DictFields(djModel,rowDict,dictList):
        for d in dictList:
            #print("fromDict",d)
            if pd.notnull(rowDict[d]):
                if d in djModel.Choice_Dictionary:
                    setattr(djModel,d,Dictionary.get(djModel.Choice_Dictionary[d],rowDict[d]))            

    choiceTables = ['COADD-Sample']
    if prgArgs.table in choiceTables:

        ExcelFile = "D:/Upload/CastDB/oraData/oraCastDB_Project_20240701_v02.xlsx"
        SheetName = "ProjectUpload"
        OutFile = "UploadProject_Issues.xlsx"

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {prgArgs.appuser}") 

        if prgArgs.table == 'Project':
            print(f"Reading {ExcelFile}.[{SheetName}]")
            df = pd.read_excel(ExcelFile,sheet_name=SheetName)
            print(df.columns)
            outDict = []
            
            # CONTACT_A_ID	CONTACT_B_ID	
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
            dictFields = ['project_type','provided_container','stock_conc_unit',]
                
            for idx,row in tqdm(df.iterrows(), total=df.shape[0]):
                
                cvPrj = Convert_ProjectID.get(row['ora_project_id'])
                if cvPrj:
                    #print(f"{row['ora_project_id']} {cvPrj.project_id}")
                    djPrj = Project.get(cvPrj.project_id)
                    if djPrj is None:
                        djPrj = Project()
                        djPrj.project_id = cvPrj.project_id
                        
                        set_fromDict(djPrj,row,cpyFields)
                        set_arrayFields_fromDict(djPrj,row,arrayFields)
                        set_DictFields(djPrj,row,dictFields)
                        if prgArgs.upload:
                            validStatus = True
                            djPrj.clean_Fields()
                            validDict = djPrj.validate()
                            if validDict:
                                validStatus = False
                                for k in validDict:
                                    print('Warning',k,validDict[k],'-')
                            if validStatus:
                                djPrj.save()
                    else:
                        row['Issue'] = f"Exists"
                        outDict.append(row)
                else:
                    row['Issue'] = f"ConvProject not found"
                    outDict.append(row)
            if len(outDict) > 0:
                print(f"Writing Issues: {OutFile}")
                outDF = pd.DataFrame(outDict)
                outDF.to_excel(OutFile)
            else:
                print(f"No Issues")
                # djPrj = Project.get(row['project_id'])
                # if djPrj is None:
                #     djPrj = Project()
                #     djPrj.project_id = row['project_id']
                #     djPrj.old_project_id = row['old_project_id']
                # if prgArgs.upload:
                #     djPrj.save()    

    
#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
