#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from zUtils import zData

import django
from django.db.models.fields.related import ManyToOneRel

#-----------------------------------------------------------------------------

import logging
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='upload_OrgDB_Data', 
                                description="Uploading data to adjCOADD from Oracle or Excel")
#    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [User]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
#    prgParser.add_argument("--excel",default=None,required=False, dest="excel", action='store', help="Excel file to upload")
#    prgParser.add_argument("-f","--vitekfolder",default=None,required=False, dest="vitekfolder", action='store', help="Vitek Folder to parse")
#    prgParser.add_argument("-p","--vitekfile",default=None,required=False, dest="vitekfile", action='store', help="Vitek File to parse")
#    prgParser.add_argument("--orgbatch",default=None,required=False, dest="orgbatch", action='store', help="OrganismBatch ID")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")
    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    djDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD"
    uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/impdata/Data"
    if prgArgs.database == 'Work':
        djDir = "I:/DEEPMICROB-Q3967/Code/Python/Django/adjCOADD"
        uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/impdata/Data"
    elif prgArgs.database == 'WorkLinux':
        djDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD"
        uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/impdata/Data"

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import ApplicationUser, Dictionary
    from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
    from apputil.utils import Validation_Log

    # Logger ----------------------------------------------------------------
    logTime= datetime.datetime.now()
    logName = "test_FK"
    logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")


    logger = logging.getLogger(__name__)
    logging.basicConfig(
#    format="%(asctime)s [%(levelname)-8s] [%(name)s] %(message)s ",
    format="[%(name)-20s] %(message)s ",
    handlers=[logging.FileHandler(logFileName,mode='w'),
              logging.StreamHandler()],
#    handlers=[logging.StreamHandler()],
    level=logging.INFO)
#    level=logging.DEBUG)

    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------

    #test_FK_OrganismID('FG_0001')

    # fg01 = Organism.get(OrgID)
    # lstFK = ['organism_batch','organims_culture']

    #print(Organism._meta.related_objects)
    #print(fg01.organism_batch_organism_id.all())

    lstModels = get_related_models(Organism)
#    lstN = get_number_related_instances(lstModels,Organism._meta.pk.name,OrgID)
    OrgID = 'FG_0001'
    lstN = get_number_related_instances(Organism,OrgID)
    print(lstN)
    # for m in lstModels:
    #     qryStr = f"{m.__name__.lower()}_organism_id"
    #     print(f"{m.__name__} : {len(m.objects.all().filter(organism_id=OrgID))}")


#-----------------------------------------------------------------------------
def get_number_related_instances(model,pk_value):
#-----------------------------------------------------------------------------
    lstModels = get_related_models(model)
    related_dict = {'Total': 0}
    for m in lstModels:
        qryDict = {model._meta.pk.name: pk_value}
        related_dict[m.__name__] = len(m.objects.all().filter(**qryDict))
        related_dict['Total'] += related_dict[m.__name__]
    return(related_dict)

#-----------------------------------------------------------------------------
def get_related_models(model):
#-----------------------------------------------------------------------------
    # print(model)
    related_models = []
    for related_object in model._meta.related_objects:
        if isinstance(related_object, ManyToOneRel):
            related_models.append(related_object.related_model)
    return related_models

#-----------------------------------------------------------------------------
def get_number_related_instances_old(related_models,pk_field,pk_value):
#-----------------------------------------------------------------------------
    related_dict = {'Total': 0}
    for m in related_models:
        qryDict = {pk_field: pk_value}
        related_dict[m.__name__] = len(m.objects.all().filter(**qryDict))
        related_dict['Total'] += related_dict[m.__name__]
    return(related_dict)


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================

