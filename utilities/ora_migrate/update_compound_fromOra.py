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

def get_oraCompound(test=0):
    from oraCastDB.oraCastDB import openCastDB

    renameCol = {
        "compound_id":      "ora_compound_id",
        "project_id":       "ora_project_id",
        "compound_type":    "ora_compound_type",
        "link_spark":       "spark_id",
        "link_chembl":      "chembl_id",
        "reg_solubility":   "reg_solvent",
    }

    replaceValues = {
      'reg_conc_unit':{'ug/ul':'mg/mL','mg/ml':'mg/mL'},
    }

    cmpSQL = "Select * From Compound Where Substr(Compound_ID,0,2) in ('C0','CX') "
    # Leaving MCC (3132), CM (190) and S00 (1) - from ora.Compound

    if test>0:
        cmpSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[Compounds] ... ")
    cmpDF = pd.DataFrame(CastDB.get_dict_list(cmpSQL))
    nTotal = len(cmpDF)
    logger.info(f"[Compounds] {nTotal} ")
    CastDB.close()

    logger.info(f"DF - Rename Columns {len(renameCol)}")
    cmpDF.rename(columns=renameCol, inplace=True)

    # logger.info(f"DF - Replace Values {len(replaceValues)}")
    # for k in replaceValues:
    #     cmpDF[k].replace(replaceValues[k],inplace=True)

    return(cmpDF)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from apputil.utils.set_data import set_arrayFields, set_dictFields, set_Dictionaries
    from dsample.models import Project, COADD_Compound, Sample
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "CompoundID" :

        cmpDF = get_oraCompound(int(prgArgs.test))
        print("--> oraCompound ---------------------------------------------------------")

        print(cmpDF.columns)
        print("-------------------------------------------------------------------------")
        OutFile = f"UpdateCompound_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"

        cpyFields = ['compound_code','compound_name','compound_description',
                     'cpoz_sn','cpoz_id','coadd_id','spark_id',
                    'reg_smiles','reg_structure','reg_mw','reg_mf',
                    'reg_amount','reg_volume','reg_conc','reg_solvent',
                    'stock_amount','prep_date',
                    'pub_status','pub_date',
                    'ora_compound_id','ora_project_id'
                    ]
        # arrayFields = {'screen_status': 'screen_status',
        #                 'report_status': 'report_status',
        #                 'compound_status': 'compound_status',
        #                 'data_status': 'data_status',
        #                 'stock_status': 'stock_status',
        #                 'pub_status': 'pub_status',
        #                 'ora_contact_ids':['contact_a_id','contact_b_id']
        #             }
        dictFields = [
                    'reg_amount_unit','reg_volume_unit','reg_conc_unit','stock_amount_unit',
                      ]


        outNumbers = {'Proc':0,'New Compounds':0,'Upload Compounds':0, 'New Samples': 0, 'Upload Samples': 0}
        outDict = []    
        for idx,row in tqdm(cmpDF.iterrows(), total=cmpDF.shape[0]):
            new_compound = False
            outNumbers['Proc'] += 1
            cvPrj = Convert_ProjectID.get(row['ora_project_id'])

            if cvPrj:

                cvCmpd = Convert_CompoundID.get(row['ora_compound_id'])
                if not cvCmpd:
                    cvCmpd = Convert_CompoundID.new_COADD_Compound_ID(row['ora_compound_id'])

                if cvCmpd:
                    #print(f"{row['ora_project_id']} {cvPrj.project_id}")
                    new_compound = False
                    djCmpd = COADD_Compound.get(cvCmpd.compound_id)
                    if djCmpd is None:
                        djCmpd = COADD_Compound()
                        djCmpd.compound_id = cvCmpd.compound_id
                        new_compound = True
                        outNumbers['New Compounds'] += 1
                    else:
                        row['Issue'] = f"Exists"

                    djPrj= Project.get(cvPrj.project_id)
                    djCmpd.project_id = djPrj
                    
                    new_sample = False
                    djSample = Sample.get(djCmpd.compound_id)
                    if djSample is None:
                        djSample = Sample()
                        djSample.sample_id = cvCmpd.compound_id
                        djSample.sample_source = 'COADD'
                        new_sample = True
                        outNumbers['New Samples'] += 1

                    set_dictFields(djCmpd,row,cpyFields)
                #     set_arrayFields(djPrj,row,arrayFields)
                    set_Dictionaries(djCmpd,row,dictFields)

                    # - Sample --------------------------------------
                    djSample.sample_code = djCmpd.compound_code

                    validStatus = True

                    djSample.clean_Fields()
                    validDict = djSample.validate()
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            print('Warning',k,validDict[k],'-')
                        outDict.append(row)

                    if validStatus:
                        if prgArgs.upload:
                            if new_sample or prgArgs.overwrite:
                                outNumbers['Upload Samples'] += 1
                                djSample.save()

                    # - Compound --------------------------------------

                    djCmpd.sample_id = djSample
                    validStatus = True

                    djCmpd.clean_Fields()
                    validDict = djCmpd.validate()
                    if validDict:
                        validStatus = False
                        for k in validDict:
                            print('Warning',k,validDict[k],'-')
                        outDict.append(row)

                    if validStatus:
                        if prgArgs.upload:
                            if new_compound or prgArgs.overwrite:
                                outNumbers['Upload Compounds'] += 1
                                djCmpd.save()
                            
                else:
                    row['Issue'] = f"ConvCompound not found"
                    print(f"[oraCompound] ConvCompound {row['ora_compound_id']} not found")
                    outDict.append(row)
            else:
                row['Issue'] = f"ConvProject not found"
                print(f"[oraCompound] ConvProject {row['ora_project_id']} not found")
                outDict.append(row)

        print(f"[oraCompound] :{outNumbers}")
        if len(outDict) > 0:
            print(f"Writing Issues: {OutFile}")
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
