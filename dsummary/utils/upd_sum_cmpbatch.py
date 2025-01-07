#
import numpy as np
import pandas as pd

from django.db.models import Q

from dsummary.models import (Summary_CmpBatch,  Summary_CmpBatch_Doseresp,  Summary_CmpBatch_Inhib,
                             Summary_Structure, Summary_Structure_Doseresp, Summary_Structure_Inhib,)
from dchem.models import Chem_Structure
from dplate.models import TestWell
from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
from ddrug.utils.bio_data import DR_Range
from adjcoadd.constants import COMPOUND_SEP

import logging
logger = logging.getLogger(__name__)

# Agg Funvtions  =======================================================================

# --------------------------------------------------------------------------------------
def get_strList(x, maxN = 10):
    if len(x) > maxN:
        _v, _c = np.unique(x, return_counts=True)
        _a = []
        for _i in range(len(_v)):
            _a.append(f"{_v[_i]} ({_c[_i]})")
        return "; ".join(_a)     
    return "; ".join(x) 

# --------------------------------------------------------------------------------------
def get_strList_unique(x):
    return ";".join(set(x)) 

# --------------------------------------------------------------------------------------
def get_nAct(x):
    return len([a for a in x if a == 'A']) 

# --------------------------------------------------------------------------------------
def get_DR_Range(x):
    return DR_Range(x)


# Summary SC Function  =======================================================================
# --------------------------------------------------------------------------------------
def pivot_sum_sc(SumType,dfSC,CmpBatchLst,StructureID,OutNumbers,
                        upload=False,overwrite=False,appuser='J.Zuegg' ):
# --------------------------------------------------------------------------------------
    OutDict = []
    CmpDict = {'sc_n_assayids' : 0, 'sc_n_actives' :0,
               'gp_n_assayids' : 0, 'gp_n_actives' :0,
               'gn_n_assayids' : 0, 'gn_n_actives' :0,
               'fg_n_assayids' : 0, 'fg_n_actives' :0,
               'cc_n_assayids' : 0, 'cc_n_actives' :0,
               'hc_n_assayids' : 0, 'hc_n_actives' :0,
               'gnm_n_assayids' : 0, 'gnm_n_actives' :0,}
    
    # Group By
    pivDF = dfSC.groupby(['assay_id']).agg({'inhibition': ['mean','max','min','std'],
                                        'mscore': ['mean','size'],
                                        'act_type': [get_strList, get_nAct ],
                                        'act_score': ['mean' ],
                                        })
    #print( pivDF.columns)
    for AssayID,row in pivDF.iterrows():
        validStatus = True
        NewEntry = False

        OutNumbers['Processed'] += 1

        if SumType == 'CmpBatch':
            djSum = Summary_CmpBatch_Inhib.get(CmpBatchLst,AssayID,Exact=True,verbose=0)
            if djSum is None:
                djSum = Summary_CmpBatch_Inhib()
                djSum.set_cmpbatch_id(CmpBatchLst)
                djSum.assay_id = AssayID
                NewEntry = True
                OutNumbers['New Entry'] += 1
        elif SumType == 'Structure':
            djSum = Summary_Structure_Inhib.get(StructureID,AssayID,verbose=0)
            if djSum is None:
                djSum = Summary_Structure_Inhib()
                djSum.structure_id = Chem_Structure.get(StructureID)
                djSum.assay_id = AssayID
                NewEntry = True
                OutNumbers['New Entry'] += 1

        djSum.act_types = row[ ('act_type','get_strList')]
        djSum.n_actives = row[ ('act_type','get_nAct')]
        djSum.n_assays = row[('mscore','size')]
        djSum.act_score_ave = row[('act_score','mean')]

        djSum.inhibition_ave = row[('inhibition','mean')]
        djSum.inhibition_std = row[('inhibition','std')]
        djSum.inhibition_min = row[('inhibition','min')]
        djSum.inhibition_max = row[('inhibition','max')]
        djSum.mscore_ave = row[('mscore','mean')]

        # n_assayids and n_actives
        CmpDict['sc_n_assayids'] += 1
        if djSum.n_actives > 0:
            CmpDict['sc_n_actives'] += 1

        if 'GP' in AssayID:
            CmpDict['gp_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['gp_n_actives'] += 1
        elif 'FG' in AssayID:
            CmpDict['fg_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['fg_n_actives'] += 1
        elif 'GN' in AssayID:
            if AssayID in Summary_CmpBatch.GNM_ASSAYS:
                CmpDict['gnm_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gnm_n_actives'] += 1
            else:
                CmpDict['gn_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gn_n_actives'] += 1
        elif 'CL' in AssayID:
            if 'CC50' in row[('result_type','get_strList_unique')]:
                CmpDict['cc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['cc_n_actives'] += 1
            if 'HC50' in row[('result_type','get_strList_unique')]:
                CmpDict['hc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['hc_n_actives'] += 1

        # Vakidate and Save
        djSum.clean_Fields()
        validDict = djSum.validate()
        if validDict:
            validStatus = False
            # for k in validDict:
            #     print('Warning',k,validDict[k],'-')
            row.update(validDict)
            OutDict.append(row)

        if validStatus:
            if upload:
                if NewEntry or overwrite:
                    #djSum.chk_migration = 0
                    OutNumbers['Upload Entries'] += 1
                    djSum.save(user=appuser)

    return(CmpDict,OutDict,OutNumbers)   

# --------------------------------------------------------------------------------------
def sum_cmpbatch_sc(CmpBatchLst,upload=False,overwrite=False, appuser='J.Zuegg'):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0,'Empty Entries':0}
    OutDict = []

    NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
    if djCmpd is None:
        djCmpd = Summary_CmpBatch()
        djCmpd.set_cmpbatch_id(CmpBatchLst)

    djCmpd.sc_n_assayids = 0
    djCmpd.sc_n_actives = 0
    djCmpd.sc_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djCmpd.sc_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)

    qryInhib = TestWell.objects.filter(cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    plate_id__result_type = 'Inhibition',
                                    is_valid = True,
                                    plate_id__plate_quality = 'Valid'
                                    ).exclude(plate_id__readout_type = 'Visual').values(
                                        'plate_id','well_id','plate_id__result_type','plate_id__assay_id',
                                        'inhibition','mscore','act_type'
                                            )

    if qryInhib.exists():
        dfInhib = pd.DataFrame(qryInhib)
        dfInhib.rename(columns={'plate_id__assay_id':'assay_id',
                             'plate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict, OutNumbers = pivot_sum_sc('CmpBatch',dfInhib,CmpBatchLst,None,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser)
        
        djCmpd.sc_n_assayids += _cmpdict['sc_n_assayids']
        djCmpd.sc_n_actives += _cmpdict['sc_n_actives']

        for _a in ['gp','gn','fg','gnm']:
            djCmpd.sc_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djCmpd.sc_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']
    else:
        OutNumbers['Empty Entries'] += 1

    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djCmpd.save(user=appuser)

    return(OutNumbers,OutDict)

# --------------------------------------------------------------------------------------
def sum_structure_sc(StructureID,upload=False,overwrite=False, appuser='J.Zuegg'):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0,'Empty Entries':0}
    OutDict = []

    #NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djStr = Summary_Structure.get(StructureID,verbose=0)
    if djStr is None:
        djStr = Summary_Structure()
        djStr.structure_id = Chem_Structure.get(StructureID)

    djStr.sc_n_assayids = 0
    djStr.sc_n_actives = 0
    djStr.sc_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djStr.sc_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)



    qryInhib = TestWell.objects.filter(cmpbatch_id__structure_id = StructureID,
                                    n_cmpbatches = 1, 
                                    plate_id__result_type = 'Inhibition',
                                    is_valid = True,
                                    plate_id__plate_quality = 'Valid'
                                    ).exclude(plate_id__readout_type = 'Visual').values(
                                        'plate_id','well_id','plate_id__result_type','plate_id__assay_id',
                                        'inhibition','mscore','act_type','act_score'
                                            )

    if qryInhib.exists():
        dfInhib = pd.DataFrame(qryInhib)
        dfInhib.rename(columns={'plate_id__assay_id':'assay_id',
                             'plate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict, OutNumbers = pivot_sum_sc('Structure',dfInhib,None,StructureID,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser)
        
        djStr.sc_n_assayids += _cmpdict['sc_n_assayids']
        djStr.sc_n_actives += _cmpdict['sc_n_actives']

        # print(f" {_cmpdict}")
        # print(f" {Summary_CmpBatch.ASSAY_CLASSES}")
        for _a in ['gp','gn','fg','gnm']:
            djStr.sc_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djStr.sc_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']
    else:
        OutNumbers['Empty Entries'] += 1
    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djStr.save(user=appuser)

    return(OutNumbers,OutDict)

# Summary DR Function  =======================================================================

# --------------------------------------------------------------------------------------
def pivot_sum_dr(SumType,drType,dfDR,CmpBatchLst,StructureID,OutNumbers,
                        upload=False,overwrite=False,appuser='J.Zuegg' ):
# --------------------------------------------------------------------------------------
    OutDict = []
    CmpDict = {'dr_n_assayids' : 0, 'dr_n_actives' :0,
               'gp_n_assayids' : 0, 'gp_n_actives' :0,
               'gn_n_assayids' : 0, 'gn_n_actives' :0,
               'fg_n_assayids' : 0, 'fg_n_actives' :0,
               'cc_n_assayids' : 0, 'cc_n_actives' :0,
               'hc_n_assayids' : 0, 'hc_n_actives' :0,
               'gnm_n_assayids' : 0, 'gnm_n_actives' :0,}

    # Group By
    pivDF = dfDR.groupby(['assay_id']).agg({'dr': [DR_Range],
                                        'inhibit_max': ['mean'],
                                        'pscore': ['mean'],
                                        'act_score': ['mean'],
                                        'act_type': [get_strList, get_nAct ],
                                        'dr_unit': [get_strList, get_strList_unique ],
                                        'result_type': [get_strList_unique ],
                                        })
    #print( pivDF.columns)
    for AssayID,row in pivDF.iterrows():
        validStatus = True
        NewEntry = False

        OutNumbers['Processed'] += 1

        if SumType == 'CmpBatch':
            djSum = Summary_CmpBatch_Doseresp.get(CmpBatchLst,AssayID,Exact=True,verbose=0)
            if djSum is None:
                djSum = Summary_CmpBatch_Doseresp()
                djSum.set_cmpbatch_id(CmpBatchLst)
                djSum.assay_id = AssayID
                NewEntry = True
        elif SumType == 'Structure':
            djSum = Summary_Structure_Doseresp.get(StructureID,AssayID,verbose=0)
            if djSum is None:
                djSum = Summary_Structure_Doseresp()
                djSum.assay_id = AssayID
                djSum.structure_id = Chem_Structure.get(StructureID)
                NewEntry = True


        djSum.drval_type = drType
        djSum.act_types = row[('act_type','get_strList')]
        djSum.n_actives = row[('act_type','get_nAct')]
        djSum.act_score_ave = row[('act_score','mean')]
        djSum.inhibit_max_ave = row[('inhibit_max','mean')]
        djSum.pscore_ave = row[('pscore','mean')]
        djSum.drval_max    = row[('dr','DR_Range')]['Max']
        djSum.drval_min    = row[('dr','DR_Range')]['Min']
        djSum.drval_median = row[('dr','DR_Range')]['Median']
        djSum.n_assays = row[('dr','DR_Range')]['nDR']
        djSum.drval_unit   = row[('dr_unit','get_strList_unique')]

        # n_assayids and n_actives
        CmpDict['dr_n_assayids'] += 1
        if djSum.n_actives > 0:
            CmpDict['dr_n_actives'] += 1

        if 'GP' in AssayID:
            CmpDict['gp_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['gp_n_actives'] += 1
        elif 'FG' in AssayID:
            CmpDict['fg_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['fg_n_actives'] += 1
        elif 'GN' in AssayID:
            if AssayID in Summary_CmpBatch.GNM_ASSAYS:
                CmpDict['gnm_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gnm_n_actives'] += 1
            else:
                CmpDict['gn_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gn_n_actives'] += 1
        elif 'CL' in AssayID:
            if 'CC50' in row[('result_type','get_strList_unique')]:
                CmpDict['cc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['cc_n_actives'] += 1
            if 'HC50' in row[('result_type','get_strList_unique')]:
                CmpDict['hc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['hc_n_actives'] += 1

        # Validate and Save
        djSum.clean_Fields()
        validDict = djSum.validate()
        if validDict:
            validStatus = False
            # for k in validDict:
            #     print('Warning',k,validDict[k],'-')
            row.update(validDict)
            #logger.warning(f"{djSum.assay_id} {djSum.structure_id} {validDict} {djSum.act_types}")
            OutDict.append(row)

        if validStatus:
            if upload:
                if NewEntry or overwrite:
                    #djSum.chk_migration = 0
                    OutNumbers['Upload Entries'] += 1
                    djSum.save(user=appuser)

    return(CmpDict,OutDict)   

# --------------------------------------------------------------------------------------
def sum_cmpbatch_dr(CmpBatchLst,upload=False,overwrite=False, appuser='J.Zuegg'):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
    if djCmpd is None:
        djCmpd = Summary_CmpBatch()
        djCmpd.set_cmpbatch_id(CmpBatchLst)

    djCmpd.dr_n_assayids = 0
    djCmpd.dr_n_actives = 0
    djCmpd.dr_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djCmpd.dr_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    # - MIC ----------------------------------------------------------
    qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryMIC.exists():
        dfDR = pd.DataFrame(qryMIC)
        dfDR.rename(columns={'testplate_id__assay_id':'assay_id','mic':'dr','mic_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('CmpBatch','MIC',dfDR,CmpBatchLst,None,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
        djCmpd.dr_n_actives += _cmpdict['dr_n_actives']

        # print(f" {_cmpdict}")
        # print(f" {Summary_CmpBatch.ASSAY_CLASSES}")
        for _a in ['gp','gn','fg','gnm']:
            djCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

   # - CC50 ----------------------------------------------------------
    qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryCC50.exists():
        dfDR = pd.DataFrame(qryCC50)
        dfDR.rename(columns={'testplate_id__assay_id':'assay_id','cc50':'dr','cc50_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('CmpBatch','CC50',dfDR,CmpBatchLst,None,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
        djCmpd.dr_n_actives += _cmpdict['dr_n_actives']
        for _a in ['cc']:
            djCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

   # - HC50 ----------------------------------------------------------
    qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','hc50','hc50_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryHC50.exists():
        dfDR = pd.DataFrame(qryHC50)
        dfDR.rename(columns={'testplate_id__assay_id':'assay_id','hc50':'dr','hc50_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('CmpBatch','HC50',dfDR,CmpBatchLst,None,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
        djCmpd.dr_n_actives += _cmpdict['dr_n_actives']
        for _a in ['hc']:
            djCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djCmpd.save(user=appuser)

    return(OutNumbers,OutDict)

# --------------------------------------------------------------------------------------
def sum_structure_dr(StructureID,upload=False,overwrite=False, appuser='J.Zuegg'):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    # CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)
    # NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djStr = Summary_Structure.get(StructureID,verbose=0)
    if djStr is None:
        djStr = Summary_Structure()
        djStr.structure_id = Chem_Structure.get(StructureID)

    djStr.dr_n_assayids = 0
    djStr.dr_n_actives = 0
    djStr.dr_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djStr.dr_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)

    # - MIC ----------------------------------------------------------
    qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_id__structure_id = StructureID, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryMIC.exists():
        dfMIC = pd.DataFrame(qryMIC)
        dfMIC.rename(columns={'testplate_id__assay_id':'assay_id','mic':'dr','mic_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('Structure','MIC',dfMIC,None,StructureID,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djStr.dr_n_assayids += _cmpdict['dr_n_assayids']
        djStr.dr_n_actives += _cmpdict['dr_n_actives']

        for _a in ['gp','gn','fg','gnm']:
            djStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

    # - CC50 ----------------------------------------------------------
    qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_id__structure_id = StructureID, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryCC50.exists():
        dfCC50 = pd.DataFrame(qryCC50)
        dfCC50.rename(columns={'testplate_id__assay_id':'assay_id','cc50':'dr','cc50_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('Structure','CC50',dfCC50,None,StructureID,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djStr.dr_n_assayids += _cmpdict['dr_n_assayids']
        djStr.dr_n_actives += _cmpdict['dr_n_actives']

        for _a in ['gp','gn','fg','gnm']:
            djStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

    # - HC50 ----------------------------------------------------------
    qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_id__structure_id = StructureID, 
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values('testplate_id__assay_id','hc50','hc50_unit','act_type','act_score','pscore','inhibit_max',
                                             'testplate_id','testwell_id','testplate_id__result_type',)

    if qryHC50.exists():
        dfHC50 = pd.DataFrame(qryHC50)
        dfHC50.rename(columns={'testplate_id__assay_id':'assay_id','hc50':'dr','hc50_unit': 'dr_unit',
                             'testplate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict = pivot_sum_dr('Structure','HC50',dfHC50,None,StructureID,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        djStr.dr_n_assayids += _cmpdict['dr_n_assayids']
        djStr.dr_n_actives += _cmpdict['dr_n_actives']

        for _a in ['gp','gn','fg','gnm']:
            djStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']

    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djStr.save(user=appuser)

    return(OutNumbers,OutDict)

    # - MIC ----------------------------------------------------------
