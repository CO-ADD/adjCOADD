#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Upload_CompoundID"
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

    cmpSQL = "Select * From Compound Where Substr(Compound_ID,0,2) in ('C0','CX','CM') And is_migrated < 1"
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
    from adjCOADD.applib.data.set_fielddata import set_model_arrayfields, set_dictFields, set_model_dicts
    from dsample.models import Project, COADD_Compound, Compound_Batch
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    # CO C0000001 - C0378034  

    # C0 C000000001 - C000400001
    # CM C000400000 - 
    # CX C010000001 - C010362007

    if prgArgs.table == "CompoundID" :

        print("--> oraCompound ---------------------------------------------------------")
        cmpDF = get_oraCompound(int(prgArgs.test))
        print(cmpDF.columns)
        print("-------------------------------------------------------------------------")
        OutFile = f"UpdateCompound_fromORA_{logTime:%Y%m%d_%H%M%S}.xlsx"

        cpyFields = ['compound_code','compound_name','compound_description',
                     'cpoz_sn','cpoz_id','coadd_id','spark_id',
                    'reg_smiles','reg_structure','reg_mw','reg_mf',
                    'reg_amount','reg_volume','reg_conc','reg_solvent',
                    'stock_amount','prep_date',
                    'pub_status','pub_date',
                    'ora_compound_id','ora_project_id','ora_compound_type'
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


        outNumbers = {'Proc':0,'New Compounds':0,'Upload Compounds':0, 'New Batches': 0, 'Upload Batches': 0}
        outDict = []    
        for idx,row in tqdm(cmpDF.iterrows(), total=cmpDF.shape[0]):
            #print(row)
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

                    # Only process new ones of to overwrite
                    if new_compound or prgArgs.overwrite:

                        djPrj= Project.get(cvPrj.project_id)
                        djCmpd.project_id = djPrj
                        
                        new_batch = False
                        djBatch = Compound_Batch.get(djCmpd.compound_id)
                        if djBatch is None:
                            djBatch = Compound_Batch()
                            djBatch.cmpbatch_id = cvCmpd.compound_id
                            djBatch.batch_source = 'COADD'
                            new_batch = True
                            outNumbers['New Batches'] += 1

                        set_dictFields(djCmpd,row,cpyFields)
                    #     set_model_arrayfields(djPrj,row,arrayFields)
                        set_model_dicts(djCmpd,row,dictFields)

                        if djCmpd.reg_mw < 2:
                            djCmpd.reg_mw = 0
                        if djCmpd.reg_mf == 'CxHxNxOx':
                            djCmpd.reg_mf = ''    

                        # - CmpBatch --------------------------------------
                        djBatch.batch_code = djCmpd.compound_code

                        validStatus = True

                        djBatch.set_defaults_model()
                        validDict = djBatch.validate_fields()
                        if validDict:
                            validStatus = False
                            for k in validDict:
                                print('Warning',k,validDict[k],'-')
                            outDict.append(row)

                        if validStatus:
                            if prgArgs.upload:
                                if new_batch or prgArgs.overwrite:
                                    outNumbers['Upload Batches'] += 1
                                    djBatch.save()

                        # - Compound --------------------------------------
                        djCmpd.cmpbatch_id = djBatch
                        validStatus = True

                        djCmpd.set_defaults_model()
                        validDict = djCmpd.validate_fields()
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
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [CompoundID]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
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
