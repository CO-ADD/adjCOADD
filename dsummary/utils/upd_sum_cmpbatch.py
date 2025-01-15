#
import numpy as np
import pandas as pd

from django.db.models import Q

from dsummary.models import (Summary_CmpBatch,  Summary_CmpBatch_Doseresp,  Summary_CmpBatch_Inhib,
                             Summary_Structure, Summary_Structure_Doseresp, Summary_Structure_Inhib,)
from dchem.models import Chem_Structure
from dplate.models import TestWell
from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run, Assay
from ddrug.utils.bio_data import DR_Range, conv_Conc, split_DR, format_DR, DR_GeoMean
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

def get_DR_GeoMean(x):
    return DR_GeoMean(x)

# --------------------------------------------------------------------------------------
def apply_DR_Std(s):
    if s['full_mw'] > 0:
        _prefix, _val, _ = split_DR(s['dr'])

        _val_uM,_unit_uM =conv_Conc(_val,s['dr_unit'],'uM',mw=s['full_mw'])
        s['dr_uM'] = format_DR(_prefix,_val_uM)
        s['dr_uM_unit'] = _unit_uM

        _val_ug,_unit_ug =conv_Conc(_val,s['dr_unit'],'ug/mL',mw=s['full_mw'])
        s['dr_ug'] = format_DR(_prefix,_val_ug)
        s['dr_ug_unit'] = _unit_ug
    else:
        s['dr_uM'] = ''
        s['dr_uM_unit'] = ''
        s['dr_ug'] = s['dr']
        s['dr_ug_unit'] = s['dr_unit']

    return(s)

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
    pivDF = dfSC.groupby(['sum_assay_id']).agg({'inhibition': ['mean','max','min','std'],
                                        'mscore': ['mean','size'],
                                        'act_type': [get_strList, get_nAct ],
                                        'act_score': ['mean' ],
                                        })
    #print( pivDF.columns)
    for SumAssayID,row in pivDF.iterrows():
        validStatus = True
        NewEntry = False

        OutNumbers['Processed'] += 1

        if SumType == 'CmpBatch':
            djSum = Summary_CmpBatch_Inhib.get(CmpBatchLst,SumAssayID,Exact=True,verbose=0)
            if djSum is None:
                djSum = Summary_CmpBatch_Inhib()
                djSum.set_cmpbatch_id(CmpBatchLst)
                djSum.sum_assay_id = SumAssayID
                NewEntry = True
                OutNumbers['New Entry'] += 1
        elif SumType == 'Structure':
            djStructure = Chem_Structure.get(StructureID)
            djSum = Summary_Structure_Inhib.get(djStructure,SumAssayID,verbose=0)
            if djSum is None:
                djSum = Summary_Structure_Inhib()
                djSum.structure_id = Chem_Structure.get(StructureID)
                djSum.sum_assay_id = SumAssayID
                NewEntry = True
                OutNumbers['New Entry'] += 1

        djSum.act_types = row[ ('act_type','get_strList')]
        djSum.n_actives = row[ ('act_type','get_nAct')]
        djSum.n_assays = row[('mscore','size')]
        #djSum.act_score_ave = row[('act_score','mean')]

        djSum.set_actscores()

        # print(" ")
        # print(f"{djSum.n_assays} {djSum.n_actives} {djSum.act_types}")

        djSum.inhibition_ave = row[('inhibition','mean')]
        if np.isnan(row[('inhibition','std')]):
            djSum.inhibition_std = 0
        else:
            djSum.inhibition_std = row[('inhibition','std')]
        djSum.inhibition_min = row[('inhibition','min')]
        djSum.inhibition_max = row[('inhibition','max')]
        djSum.mscore_ave = row[('mscore','mean')]

        # n_assayids and n_actives
        CmpDict['sc_n_assayids'] += 1
        if djSum.n_actives > 0:
            CmpDict['sc_n_actives'] += 1

        if 'GP' in SumAssayID:
            CmpDict['gp_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['gp_n_actives'] += 1
        elif 'FG' in SumAssayID:
            CmpDict['fg_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['fg_n_actives'] += 1
        elif 'GN' in SumAssayID:
            if SumAssayID in Summary_CmpBatch.GNM_ASSAYS:
                CmpDict['gnm_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gnm_n_actives'] += 1
            else:
                CmpDict['gn_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gn_n_actives'] += 1
        elif 'CL' in SumAssayID:
            if 'CC50' in row[('result_type','get_strList_unique')]:
                CmpDict['cc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['cc_n_actives'] += 1
            if 'HC50' in row[('result_type','get_strList_unique')]:
                CmpDict['hc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['hc_n_actives'] += 1

        # Vakidate and Save
        djSum.init_fields()
        validDict = djSum.validate_fields()
        if validDict:
            validStatus = False
            # for k in validDict:
            #     print('Warning',k,validDict[k],'-')
            row.update(validDict)
            logger.warning(f" {validDict}")
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
    if NCmpBatches > 0:
        # Sum_Cmpd ---------------------------------------------------------------
        djSumCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
        if djSumCmpd is None:
            djSumCmpd = Summary_CmpBatch()
            djSumCmpd.set_cmpbatch_id(CmpBatchLst)

        djSumCmpd.sc_n_assayids = 0
        djSumCmpd.sc_n_actives = 0
        djSumCmpd.sc_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
        djSumCmpd.sc_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)

        qryInhib = TestWell.objects.filter(cmpbatch_lst__contains = CmpBatchLst, 
                                        n_cmpbatches = NCmpBatches, 
                                        plate_id__result_type = 'Inhibition',
                                        is_valid = True,
                                        plate_id__plate_quality = 'Valid'
                                        ).exclude(plate_id__readout_type = 'Visual').values(
                                            'plate_id','well_id','plate_id__result_type','plate_id__assay_id__sum_assay_id',
                                            'inhibition','mscore','act_type','act_score'
                                                )

        if qryInhib.exists():
            dfInhib = pd.DataFrame(qryInhib)
            dfInhib.rename(columns={'plate_id__assay_id__sum_assay_id':'sum_assay_id',
                                'plate_id__result_type':'result_type',}, inplace=True)

            _cmpdict, _outdict, OutNumbers = pivot_sum_sc('CmpBatch',dfInhib,CmpBatchLst,None,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser)
            
            djSumCmpd.sc_n_assayids += _cmpdict['sc_n_assayids']
            djSumCmpd.sc_n_actives += _cmpdict['sc_n_actives']

            for _a in ['gp','gn','fg','gnm']:
                djSumCmpd.sc_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
                djSumCmpd.sc_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']
        else:
            OutNumbers['Empty Entries'] += 1

        # Sum_Cmpd ---------------------------------------------------------------
        if upload:
            djSumCmpd.save(user=appuser)
    else:
        OutNumbers['Empty Entries'] += 1
    return(OutNumbers,OutDict)

# --------------------------------------------------------------------------------------
def sum_structure_sc(StructureID,upload=False,overwrite=False, appuser='J.Zuegg'):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0,'Empty Entries':0}
    OutDict = []

    #NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djSumStr = Summary_Structure.get(StructureID,verbose=0)
    if djSumStr is None:
        djSumStr = Summary_Structure()
        djSumStr.structure_id = Chem_Structure.get(StructureID)

    djSumStr.sc_n_assayids = 0
    djSumStr.sc_n_actives = 0
    djSumStr.sc_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djSumStr.sc_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)



    qryInhib = TestWell.objects.filter(cmpbatch_id__structure_id = StructureID,
                                    n_cmpbatches = 1, 
                                    plate_id__result_type = 'Inhibition',
                                    is_valid = True,
                                    plate_id__plate_quality = 'Valid'
                                    ).exclude(plate_id__readout_type = 'Visual').values(
                                        'plate_id','well_id','plate_id__result_type','plate_id__assay_id__sum_assay_id',
                                        'inhibition','mscore','act_type','act_score'
                                            )

    if qryInhib.exists():
        dfInhib = pd.DataFrame(qryInhib)
        dfInhib.rename(columns={'plate_id__assay_id__sum_assay_id':'sum_assay_id',
                             'plate_id__result_type':'result_type',}, inplace=True)
        _cmpdict, _outdict, OutNumbers = pivot_sum_sc('Structure',dfInhib,None,StructureID,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser)
        
        djSumStr.sc_n_assayids += _cmpdict['sc_n_assayids']
        djSumStr.sc_n_actives += _cmpdict['sc_n_actives']

        # print(f" {_cmpdict}")
        # print(f" {Summary_CmpBatch.ASSAY_CLASSES}")
        for _a in ['gp','gn','fg','gnm']:
            djSumStr.sc_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_assayids']
            djSumStr.sc_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] = _cmpdict[f'{_a}_n_actives']
    else:
        OutNumbers['Empty Entries'] += 1
    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djSumStr.save(user=appuser)

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
    
    # Add dr_std [uM] - dr_ug, dr_ug_unit, dr_uM, dr_uM_unit
    dfDR = dfDR.apply(apply_DR_Std,axis=1)

    # Group By
    pivDF = dfDR.groupby(['sum_assay_id']).agg({'inhibit_max': ['mean'],
                                        'dr_ug': [DR_Range],
                                        'dr_uM': [DR_GeoMean],
                                        'dr_unit' : [get_strList_unique],
                                        'dr_ug_unit': [get_strList_unique ],
                                        'dr_uM_unit': [get_strList_unique ],        
                                        'pscore': ['mean'],
                                        'act_score': ['mean'],
                                        'act_type': [get_strList, get_nAct ],
                                        'result_type': [get_strList_unique ],
                                        })
    #print( pivDF.columns)
    for SumAssayID,row in pivDF.iterrows():
        validStatus = True
        NewEntry = False

        OutNumbers['Processed'] += 1
        if SumType == 'CmpBatch':            
            djSum = Summary_CmpBatch_Doseresp.get(CmpBatchLst,SumAssayID,Exact=True,verbose=0)
            if djSum is None:
                djSum = Summary_CmpBatch_Doseresp()
                djSum.set_cmpbatch_id(CmpBatchLst)
                djSum.sum_assay_id = SumAssayID
                NewEntry = True
        elif SumType == 'Structure':
            djStructure = Chem_Structure.get(StructureID)
            djSum = Summary_Structure_Doseresp.get(djStructure,SumAssayID,verbose=0)
            if djSum is None:
                djSum = Summary_Structure_Doseresp()
                djSum.sum_assay_id = SumAssayID
                djSum.structure_id = djStructure
                NewEntry = True


        djSum.drval_type = drType
        djSum.act_types = row[('act_type','get_strList')]
        djSum.n_actives = row[('act_type','get_nAct')]
        djSum.inhibit_max_ave = row[('inhibit_max','mean')]

        djSum.drval_max    = row[('dr_ug','DR_Range')]['Max']
        djSum.drval_min    = row[('dr_ug','DR_Range')]['Min']
        djSum.drval_median = row[('dr_ug','DR_Range')]['Median']
        djSum.drval_unit   = row[('dr_ug_unit','get_strList_unique')]
        djSum.n_assays = row[('dr_ug','DR_Range')]['nDR']

        djSum.drval_std_unit   = row[('dr_uM_unit','get_strList_unique')]
        djSum.drval_std_geomean = format_DR(split_DR(djSum.drval_median)[0],
                                            row[('dr_uM','DR_GeoMean')])
        djSum.set_actscores()

        #print(f" {djSum.structure_id} {djSum.drval_std_geomean} {djSum.drval_std_unit} ")
        # n_assayids and n_actives
        CmpDict['dr_n_assayids'] += 1
        if djSum.n_actives > 0:
            CmpDict['dr_n_actives'] += 1
        
        if 'GP' in SumAssayID:
            CmpDict['gp_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['gp_n_actives'] += 1
        elif 'FG' in SumAssayID:
            CmpDict['fg_n_assayids'] += 1
            if djSum.n_actives > 0:
                CmpDict['fg_n_actives'] += 1
        elif 'GN' in SumAssayID:
            if SumAssayID in Summary_CmpBatch.GNM_ASSAYS:
                CmpDict['gnm_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gnm_n_actives'] += 1
            else:
                CmpDict['gn_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['gn_n_actives'] += 1
        elif 'CL' in SumAssayID:
            if 'CC50' in row[('result_type','get_strList_unique')]:
                CmpDict['cc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['cc_n_actives'] += 1
            if 'HC50' in row[('result_type','get_strList_unique')]:
                CmpDict['hc_n_assayids'] += 1
                if djSum.n_actives > 0:
                    CmpDict['hc_n_actives'] += 1

        # Validate and Save
        djSum.init_fields()
        validDict = djSum.validate_fields()
        if validDict:
            validStatus = False
            # for k in validDict:
            #     print('Warning',k,validDict[k],'-')
            row.update(validDict)
            #logger.warning(f"{djSum.sum_assay_id} {djSum.structure_id} {validDict} {djSum.act_types}")
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
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0, 'Empty Entries':0}
    OutDict = []

    NCmpBatches = len(CmpBatchLst)

    if NCmpBatches > 0:
        # Sum_Cmpd ---------------------------------------------------------------
        djSumCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
        if djSumCmpd is None:
            djSumCmpd = Summary_CmpBatch()
            djSumCmpd.set_cmpbatch_id(CmpBatchLst)

        djSumCmpd.dr_n_assayids = 0
        djSumCmpd.dr_n_actives = 0
        djSumCmpd.dr_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
        djSumCmpd.dr_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)

        # - MIC ----------------------------------------------------------
        qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_lst__contains = CmpBatchLst, 
                                        n_cmpbatches = NCmpBatches, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')
        EmptyEntry = True
        if qryMIC.exists():
            EmptyEntry = False
            dfDR = pd.DataFrame(qryMIC)

            dfDR.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','mic':'dr','mic_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)

            _cmpdict, _outdict = pivot_sum_dr('CmpBatch','MIC',dfDR,CmpBatchLst,None,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumCmpd.dr_n_actives += _cmpdict['dr_n_actives']

            # print(f" {_cmpdict}")
            # print(f" {Summary_CmpBatch.ASSAY_CLASSES}")
            for _a in ['gp','gn','fg','gnm']:
                djSumCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

        # - CC50 ----------------------------------------------------------
        qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_lst__contains = CmpBatchLst, 
                                        n_cmpbatches = NCmpBatches, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')

        if qryCC50.exists():
            EmptyEntry = False
            dfDR = pd.DataFrame(qryCC50)
            dfDR.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','cc50':'dr','cc50_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)
            _cmpdict, _outdict = pivot_sum_dr('CmpBatch','CC50',dfDR,CmpBatchLst,None,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumCmpd.dr_n_actives += _cmpdict['dr_n_actives']
            for _a in ['cc']:
                djSumCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

        # - HC50 ----------------------------------------------------------
        qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_lst__contains = CmpBatchLst, 
                                        n_cmpbatches = NCmpBatches, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','hc50','hc50_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')

        if qryHC50.exists():
            EmptyEntry = False
            dfDR = pd.DataFrame(qryHC50)
            dfDR.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','hc50':'dr','hc50_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)
            _cmpdict, _outdict = pivot_sum_dr('CmpBatch','HC50',dfDR,CmpBatchLst,None,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumCmpd.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumCmpd.dr_n_actives += _cmpdict['dr_n_actives']
            for _a in ['hc']:
                djSumCmpd.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumCmpd.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

        # Sum_Cmpd ---------------------------------------------------------------
        if upload:
            djSumCmpd.save(user=appuser)
        
        if EmptyEntry:
            OutNumbers['Empty Entries'] += 1

    else:
        OutNumbers['Empty Entries'] += 1
    return(OutNumbers,OutDict)

# --------------------------------------------------------------------------------------
def sum_structure_dr(StructureID,upload=False,overwrite=False, appuser='J.Zuegg',AssayData=['MIC','CC50','HC50']):
# --------------------------------------------------------------------------------------
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    # CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)
    # NCmpBatches = len(CmpBatchLst)


    # Sum_Structure ---------------------------------------------------------------
    djSumStr = Summary_Structure.get(StructureID,verbose=0)
    if djSumStr is None:
        djSumStr = Summary_Structure()
        djSumStr.structure_id = Chem_Structure.get(StructureID)

    djSumStr.dr_n_assayids = 0
    djSumStr.dr_n_actives = 0
    djSumStr.dr_assayid_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)
    djSumStr.dr_actives_lst = [0] * len(Summary_CmpBatch.ASSAY_CLASSES)

    # - MIC ----------------------------------------------------------
    if 'MIC' in AssayData:
        qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_id__structure_id = StructureID, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')

        if qryMIC.exists():
            dfMIC = pd.DataFrame(qryMIC)
            dfMIC.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','mic':'dr','mic_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)
            
            _cmpdict, _outdict = pivot_sum_dr('Structure','MIC',dfMIC,None,StructureID,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumStr.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumStr.dr_n_actives += _cmpdict['dr_n_actives']
            for _a in ['gp','gn','fg','gnm']:
                djSumStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

    # - CC50 ----------------------------------------------------------
    if 'CC50' in AssayData:
        qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_id__structure_id = StructureID, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')

        if qryCC50.exists():
            dfCC50 = pd.DataFrame(qryCC50)
            dfCC50.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','cc50':'dr','cc50_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)
            _cmpdict, _outdict = pivot_sum_dr('Structure','CC50',dfCC50,None,StructureID,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumStr.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumStr.dr_n_actives += _cmpdict['dr_n_actives']

            for _a in ['gp','gn','fg','gnm']:
                djSumStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

    # - HC50 ----------------------------------------------------------
    if 'HC50' in AssayData:
        qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                        cmpbatch_id__structure_id = StructureID, 
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values('testplate_id__assay_id__sum_assay_id','hc50','hc50_unit','act_type','act_score','pscore','inhibit_max',
                                                'testplate_id','testwell_id','testplate_id__result_type',
                                                'cmpbatch_id__full_mw')

        if qryHC50.exists():
            dfHC50 = pd.DataFrame(qryHC50)
            dfHC50.rename(columns={'testplate_id__assay_id__sum_assay_id':'sum_assay_id','hc50':'dr','hc50_unit': 'dr_unit',
                                'testplate_id__result_type':'result_type','cmpbatch_id__full_mw':'full_mw'}, inplace=True)
            _cmpdict, _outdict = pivot_sum_dr('Structure','HC50',dfHC50,None,StructureID,OutNumbers,
                                                upload=upload,overwrite=overwrite,appuser=appuser )
            
            djSumStr.dr_n_assayids += _cmpdict['dr_n_assayids']
            djSumStr.dr_n_actives += _cmpdict['dr_n_actives']

            for _a in ['gp','gn','fg','gnm']:
                djSumStr.dr_assayid_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_assayids']
                djSumStr.dr_actives_lst[Summary_CmpBatch.ASSAY_CLASSES[_a]] += _cmpdict[f'{_a}_n_actives']

    # Sum_Cmpd ---------------------------------------------------------------
    if upload:
        djSumStr.save(user=appuser)

    return(OutNumbers,OutDict)

    # - MIC ----------------------------------------------------------
