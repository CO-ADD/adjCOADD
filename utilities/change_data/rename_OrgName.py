#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

# from zUtils import zData

import django
#from djCOADD import djOrgDB
# from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

import logging
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main():
    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='change_OrgDB_Data', 
                                description="Changing data to adjCOADD from Excel")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    # prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
    # prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
    # prgParser.add_argument("-o","--orgbatch",default=None,required=False, dest="orgbatch", action='store', help="OrganismBatch ID")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
    # prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD"
    uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"

    xlFiles = {
        'Application': "ApplicationData_v05.xlsx",
        'Drug': "DrugData_v04.xlsx",
        'MIC': "LMIC_Data_v06.xlsx",
        'OrgDB': "OrgDB_v20_30Jun2023.xlsx",
    }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser
    import b_change_dOrganism as dOrg

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "ChangeOrgDB"
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

 
 
    choiceTables = ['Rename_OrgName','Rename_OrgID','Rename_OrgBatchID','Manual'
                    ]
    if prgArgs.table in choiceTables:

        appuser = ApplicationUser.get(prgArgs.appuser)

        logger.info(f"[Upd_djCOADD] Table: {prgArgs.table}") 
        logger.info(f"[Upd_djCOADD] User:  {appuser}") 

        if prgArgs.table == 'Rename_OrgName' and prgArgs.file:
            dOrg.rename_OrgName_xls(prgArgs.file,XlsSheet="New OrgName",
                                    upload=prgArgs.upload,uploaduser=appuser)

        if prgArgs.table == 'Rename_OrgID' and prgArgs.file:
            dOrg.rename_OrgID_xls(prgArgs.file,XlsSheet="New GP OrgID",
                                    upload=prgArgs.upload,uploaduser=appuser)
        
        if prgArgs.table == 'Manual':
            # dOrg.rename_OrgID_allBatches('GN_0696','Bacillus velezensis',None,{'strain_code':'Ba.ve BRf 79136'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_0697','Bacillus safensis',None,{'strain_code':'Ba.sa BRf 79197'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1034','Enterococcus faecium',None,{'strain_code':'En.fam PK 38'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1044','Enterococcus faecalis',None,{'strain_code':'En.fas PK 48'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1045','Enterococcus faecalis',None,{'strain_code':'En.fas PK 49'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1052','Enterococcus faecalis',None,{'strain_code':'En.fas PK 56'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1060','Enterococcus faecalis',None,{'strain_code':'En.fas PK 65'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1069','Enterococcus faecalis',None,{'strain_code':'En.fas PK 74'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1083','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1084','Enterococcus faecium',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1085','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1086','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1087','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1092','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1096','Enterococcus faecalis',None,{'strain_code':'En.fas NP U8546'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1099','Enterococcus faecium',None,{'strain_code':'En.fam NP U292'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1100','Enterococcus faecium',None,{'strain_code':'En.fam NP U299'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1109','Enterococcus faecalis',None,{'strain_code':'En.fas NP U1057'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1110','Staphylococcus aureus',None,{'strain_code':'Sa NP P5605'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1111','Enterococcus faecalis',None,{'strain_code':'En.fas NP P3721'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1113','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1120','Enterococcus faecalis',None,{'strain_code':'En.fas NP P5089'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1122','Enterococcus faecium',None,{'strain_code':'En.fam NP P4963'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1131','Enterococcus faecalis',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1134','Staphylococcus epidermidis',None,{'strain_code':'St.ep NP P4992'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1204','Staphylococcus gallinarum',None,{'strain_code':'St.ga NG 761830041'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1178','Bacillus sp.',None,{'strain_code':'Ba.sp NG 761831767'},uploaduser=appuser)
            # dOrg.rename_OrgID_allBatches('GN_1192','Staphylococcus haemolyticus',None,{'strain_code':'St.ha NG 761832918'},uploaduser=appuser)


            # dOrg.rename_OrgID_sngBatches('GN_0952_02','Acinetobacter baumannii',None,{'strain_code':'Ab IHMA 2007771'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_0981_01','Serratia marcescens',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_0982_02','Serratia marcescens',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1149_02','Escherichia coli',None,{'strain_code':'Ec EG 6'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1254_02','Escherichia coli',None,{'strain_code':'Ec ZM 2911NK-1'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1258_01','Klebsiella pneumoniae',None,{'strain_code':'Kp ZM 2716BN-1'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1270_03','Escherichia coli',None,{'strain_code':'Ec ZM 2315NKU'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1274_02','Acinetobacter baumannii',None,{'strain_code':'Ab ZM 217NK-1'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1275_02','Escherichia coli',None,{'strain_code':'Ec ZM 2118N-1'},uploaduser=appuser)

            # dOrg.rename_OrgID_sngBatches('GN_1282_02','Enterobacter cloacae',None,{'strain_code':'En.cl ZM 2491BL83A'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1285_02','Acinetobacter junii',None,{'strain_code':'Ac.ju ZM 2917BL93'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1308_01','Achromobacter xylosoxidans',None,{'strain_code':''},uploaduser=appuser)

            # dOrg.rename_OrgID_sngBatches('GN_1307_01','Pseudomonas putida',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1316_02','Pseudomonas putida',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1318_01','Pseudomonas aeruginosa',None,{'strain_code':''},uploaduser=appuser)

            # dOrg.rename_OrgID_sngBatches('GN_1323_02','Pseudomonas putida',None,{'strain_code':''},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1334_02','Pseudomonas aeruginosa',None,{'strain_code':''},uploaduser=appuser)

            # dOrg.rename_OrgID_sngBatches('GN_1115_03','Enterococcus faecalis',None,{'strain_code':'En.fas NP P5010'},uploaduser=appuser)

            # dOrg.rename_OrgID_sngBatches('GN_1207_02','Staphylococcus gallinarum',None,{'strain_code':'St.ga NG 761831834'},uploaduser=appuser)
            # dOrg.rename_OrgID_sngBatches('GN_1270_02','Enterococcus faecalis',None,{'strain_code':'En.fas ZM 2315NKU'},uploaduser=appuser)


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================
