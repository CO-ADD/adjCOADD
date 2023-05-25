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

import logging
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
def main():

    # ArgParser -------------------------------------------------------------
    prgParser = argparse.ArgumentParser(prog='test_breakpoint', description="Test for BreakPoint")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
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

    xlFiles = {
        'Application': "ApplicationData_v03.xlsx",
        'Drug': "DrugData_v02.xlsx",
    }

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    import ddrug.utils.bio_analysis as ba
    print("...")
    objBP = ba.get_BreakPoint("Tigecycline","Escherichia coli","MIC")
    print(repr(objBP))
    #print(df)

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")

#==============================================================================