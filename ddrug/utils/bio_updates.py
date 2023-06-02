import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

#import django

from ddrug.models import Drug, MIC_COADD, MIC_Pub, Breakpoint
from ddrug.utils.bio_analysis import calc_Breakpoint

from apputil.utils.data import Timer

#-------------------------------------------------------------------------------------------
def update_Breakpoint_Profile(table='MIC_COADD',upload=False, overwrite=False, OutputN=500):
    """
    Updates bp_profile and bp_source for Antibiogram (MIC_COADD) or Public data (MIC_PUB)
    """
#-------------------------------------------------------------------------------------------

    if table == 'MIC_COADD':
        if overwrite:
            djMICLst = MIC_COADD.objects.all()
        else:
            djMICLst = MIC_COADD.objects.filter(bp_profile="")
        nTotal = len(djMICLst)

        nTime  = Timer(nTotal)
        nProcessed = 0
        print(f"[    Start/ {nTotal:8d}]  ")

        for djMIC in djMICLst:

            _drugname = djMIC.drug_id.drug_name
            _orgname = djMIC.orgbatch_id.organism_id.organism_name.organism_name
            BP = calc_Breakpoint(_drugname,_orgname,"MIC",djMIC.mic)

            if upload:
                djMIC.bp_profile = BP[0]
                djMIC.bp_source = BP[1]
                djMIC.save()

            nProcessed = nProcessed + 1
            if nProcessed%OutputN == 0:
                eTime,sTime = nTime.remains(nProcessed)
                print(f"[{nProcessed:8d} / {nTotal:8d}] {eTime} >r {djMIC} {BP}")

