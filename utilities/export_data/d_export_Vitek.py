import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

import django
#from zUtils import zData

from apputil.models import ApplicationUser, Dictionary, ApplicationLog
from apputil.utils import validation_log
from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture, OrgBatch_Stock
from ddrug.models import Drug, VITEK_Card, VITEK_ID, VITEK_AST
from ddrug.utils.import_drug import *

import ddrug.utils.vitek as Vitek

def agg_ListStr(x):
    sx = ', '.join(['{:s}'.format(i) for i in x if i is not None])
    return(sx)

#-----------------------------------------------------------------------------------
def export_Vitek(OutDir):
#-----------------------------------------------------------------------------------

    AuditFields = ['astatus','acreated_at','acreated_id','aupdated_at','aupdated_id','adeleted_at','adeleted_id']

    vIDs = VITEK_ID.objects.filter(astatus__gte=0).values().annotate(
        batch_id = F("card_barcode_id__orgbatch_id__batch_id"),
        organism_id = F("card_barcode_id__orgbatch_id__organism_id"),
        orgname = F("card_barcode__orgbatch_id__organism_id__organism_name"),
        card_type = F("card_barcode__card_type"),
        card_code = F("card_barcode__card_code"),
        )
    logger.info(f"[Vitek    ] {len(vIDs)} ID's ")

    for v in vIDs:
        v['piv_value'] = f"{v['id_organism']} ({v['id_probability']})"

    dfID = pd.DataFrame(vIDs).drop(columns=AuditFields)
    pivID = dfID.pivot_table(index=['organism_id','orgname','batch_id'],columns=['card_code'],values=['piv_value'],aggfunc=agg_ListStr)

    xlDir = os.path.join(OutDir,'Vitek')
    if not os.path.exists(xlDir):
        os.makedirs(xlDir)
    xlFile = os.path.join(xlDir,"Vitek_Data.xlsx")
    with pd.ExcelWriter(xlFile) as writer:
        dfID.to_excel(writer, sheet_name='VitekID')
        pivID.to_excel(writer, sheet_name='ID')

    #print(fIDs)
