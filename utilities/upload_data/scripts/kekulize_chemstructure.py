import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import argparse

from tqdm import tqdm

import django

def main():
    from dchem.models import Chem_Structure
    from rdkit import Chem

    AllStr = Chem_Structure.objects.all()
    for s in tqdm(AllStr):        
        doSave = False
        try:
            Chem.Kekulize(s.smol)
            doSave = True
        except:
            doSave = False
            print(f"{s.structure_id}: Kekule Failed")
        if doSave:
            s.save()
    print("Hey")

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    
    # Django -------------------------------------------------------------
    djDir = "D:/Code/zdjCode/adjCOADD"
    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    # Logger ----------------------------------------------------------------
    import logging
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

    main()

    print("-------------------------------------------------------------------")

#==============================================================================
