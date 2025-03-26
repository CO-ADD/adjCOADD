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
from django.apps import apps
from django.db.models import ForeignKey, Model

#from djCOADD import djOrgDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Check_ForeignKeys"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)    
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
def get_Fields_byForeignKey(fkModel, fkField=None):
#-----------------------------------------------------------------------------------
    all_models = apps.get_models()
    
    objModel = None
    if isinstance(fkModel,Model):
        objModel = Model
    if isinstance(fkModel,str):
        for m in all_models:
            if m.__name__ == fkModel:
                objModel = m
    
    # Find models with Author as a foreign key
    models_with_fk = []
    for m in all_models:
        #print(m.__name__)
        for f in m._meta.fields:
            #print(type(f), f.related_model)
            if isinstance(f, ForeignKey) and f.related_model == objModel:
                if fkField:
                    if fkField in f.name:
                        models_with_fk.append(f)
                else:
                    models_with_fk.append(f)    

    return(models_with_fk)


#-----------------------------------------------------------------------------------
def main(prgArgs,djDir):
#-----------------------------------------------------------------------------------

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()


    from apputil.models import Dictionary
    from adjCOADD.applib.data.set_fielddata import set_model_arrayfields, set_dictFields, set_model_dicts
    from dsample.models import Project
    from dsample.models import Convert_ProjectID, Convert_CompoundID

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------
    if prgArgs.model : 
        # Find all Models with FK to OrgBatch
        fk_models = get_Fields_byForeignKey(prgArgs.model,prgArgs.field)
        print(f" {prgArgs.model} - {prgArgs.field}")
        for fk in fk_models:
            print(f" {prgArgs.model} <- {fk}")

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-m","--model",default=None,required=True, dest="model", action='store', help="ForeignKey Model to check")
    prgParser.add_argument("-f","--field",default=None,required=True, dest="field", action='store', help="ForeignKey Field to check")

#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")

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
