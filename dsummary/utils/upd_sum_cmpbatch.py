#
import numpy as np
import pandas as pd

from django.db.models import Q

from dsummary.models import Summary_CmpBatch, Summary_CmpBatch_Doseresp, Summary_CmpBatch_Inhib
from dplate.models import TestWell
from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run
from ddrug.utils.bio_data import DR_Range
from adjcoadd.constants import COMPOUND_SEP


# Agg Funvtions  -------------------------------------------------------------
def get_strList(x, maxN = 10):
    if len(x) > maxN:
        _v, _c = np.unique(x, return_counts=True)
        _a = []
        for _i in range(len(_v)):
            _a.append(f"{_v[_i]} ({_c[_i]})")
        return "; ".join(_a)     
    return "; ".join(x) 

def get_strList_unique(x):
    return ";".join(set(x)) 

def get_nAct(x):
    return len([a for a in x if a == 'A']) 

def get_DR_Range(x):
    return DR_Range(x)



# Summary DR Function  =======================================================================
def _process_sum_cmpbatch_dr(drType,qryMIC,CmpBatchLst,NCmpBatches,CmpBatchs,djCmpd,OutNumbers,
                            upload=False,overwrite=False,appuser='J.Zuegg' ):
# =========================================================================================
    OutDict = []

    dfDR = pd.DataFrame(qryMIC).assign(cmpbatchs=CmpBatchs)
    dfDR.columns = ['plate_id','well_id','result_type','assay_id',
                    'dr','dr_unit','act_type','act_score','pscore',
                    'inhibit_max',
                    'cmpbatchs']

    pivDF = dfDR.groupby(['assay_id']).agg({'dr': [DR_Range],
                                        'inhibit_max': ['mean'],
                                        'pscore': ['mean'],
                                        'act_score': ['mean'],
                                        'act_type': [get_strList, get_nAct ],
                                        'dr_unit': [get_strList, get_strList_unique ],
                                        })
    #print( pivDF.columns)
    for AssayID,row in pivDF.iterrows():
        validStatus = True
        NewEntry = False

        OutNumbers['Processed'] += 1
        djCmpd.dr_n_assayids += 1

        djSum = Summary_CmpBatch_Doseresp.get(CmpBatchLst,AssayID,Exact=True,verbose=0)
        if djSum is None:
            djSum = Summary_CmpBatch_Doseresp()
            djSum.cmpbatch_lst = CmpBatchLst
            djSum.n_cmpbatches = NCmpBatches
            djSum.assay_id = AssayID
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

        if djSum.n_actives > 0:
            djCmpd.dr_n_actives += 1

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

# Summary DR Function  =======================================================================
def sum_cmpbatch_doseresponse(CmpBatchLst,upload=False,overwrite=False,appuser='J.Zuegg'):
# =========================================================================================
    validStatus = True
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)
    NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
    if djCmpd is None:
        djCmpd = Summary_CmpBatch()
        djCmpd.cmpbatch_lst = CmpBatchLst
        djCmpd.n_cmpbatches = len(CmpBatchLst)

    djCmpd.sc_n_assayids = 0
    djCmpd.sc_n_actives = 0

    qryTW = TestWell.objects.filter(cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    plate_id__result_type = 'Inhibition',
                                    is_valid = True,
                                    plate_id__plate_quality = 'Valid'
                                    ).values(
                                        'plate_id','well_id','plate_id__result_type','plate_id__assay_id',
                                        'inhibition','mscore','act_type'
                                            )
    if qryTW.exists():
        dfSC = pd.DataFrame(qryTW).assign(cmpbatchs=CmpBatchs)
        dfSC.columns = ['plate_id','well_id','result_type','assay_id','inhibition','mscore','act_type','cmpbatchs']

        pivDF = dfSC.groupby(['assay_id']).agg({'inhibition': ['mean','max','min','std'],
                                            'mscore': ['mean','size'],
                                            'act_type': [get_strList, get_nAct ],
                                            })
        
        #print( pivDF.columns)
        for idx,row in pivDF.iterrows():
            OutNumbers['Processed'] += 1
            djCmpd.sc_n_assayids += 1

            NewEntry = False
            djSum = Summary_CmpBatch_Inhib.get(CmpBatchLst,idx,Exact=True,verbose=0)
            if djSum is None:
                djSum = Summary_CmpBatch_Inhib()
                djSum.cmpbatch_lst = CmpBatchLst
                djSum.n_cmpbatches = NCmpBatches
                djSum.assay_id = idx
                NewEntry = True
            djSum.act_types = row[ ('act_type','get_strList')]

            djSum.n_assays = row[('mscore','size')]
            djSum.n_actives = row[ ('act_type','get_nAct')]
            #djSum.act_score_ave =

            djSum.inhibition_ave = row[('inhibition','mean')]
            djSum.inhibition_std = row[('inhibition','std')]
            djSum.inhibition_min = row[('inhibition','min')]
            djSum.inhibition_max = row[('inhibition','max')]
            djSum.mscore_ave = row[('mscore','mean')]

            if djSum.n_actives > 0:
                djCmpd.sc_n_actives += 1

            djSum.clean_Fields()
            validDict = djSum.validate()
            if validDict:
                validStatus = False
                # for k in validDict:
                #     print('Warning',k,validDict[k],'-')
                row.update(validDict)

            if validStatus:
                if upload:
                    if NewEntry or overwrite:
                        #djSum.chk_migration = 0
                        OutNumbers['Upload Entries'] += 1
                        djSum.save(user=appuser)
    # Sum_Cmpd ---------------------------------------------------------------
    # djCmpd.sc_assayid_lst  = 
    # djCmpd.sc_actives_lst  = 
    if upload:
        djCmpd.save()

    return(OutNumbers,OutDict)

# Summary DR from  AssayData MIC  ========================================================
def sum_cmpbatch_doseresponse_MIC(CmpBatchLst,NCmpBatches,djCmpd,
                                  upload=False,overwrite=False,appuser='J.Zuegg'):
# =========================================================================================
    validStatus = True
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)

    qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__result_type = 'MIC',
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values(
                                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
                                        'mic','mic_unit','act_type','act_score','pscore',
                                        'inhibit_max'
                                            )
    if qryMIC.exists():
        dfDR = pd.DataFrame(qryMIC).assign(cmpbatchs=CmpBatchs)
        dfDR.columns = ['plate_id','well_id','result_type','assay_id',
                        'dr','dr_unit','act_type','act_score','pscore',
                        'inhibit_max',
                        'cmpbatchs']

        OutDict = _process_sum_cmpbatch_dr('MIC',qryMIC,CmpBatchLst,NCmpBatches,CmpBatchs,djCmpd,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )
        
        # pivDF = dfDR.groupby(['assay_id']).agg({'dr': [DR_Range],
        #                                     'inhibit_max': ['mean'],
        #                                     'pscore': ['mean'],
        #                                     'act_score': ['mean'],
        #                                     'act_type': [get_strList, get_nAct ],
        #                                     'dr_unit': [get_strList, get_strList_unique ],
        #                                     })
        # #print( pivDF.columns)
        # for AssayID,row in pivDF.iterrows():
        #     OutNumbers['Processed'] += 1
        #     djCmpd.dr_n_assayids += 1

        #     NewEntry = False
        #     djSum = Summary_CmpBatch_Doseresp.get(CmpBatchLst,AssayID,Exact=True,verbose=0)
        #     if djSum is None:
        #         djSum = Summary_CmpBatch_Doseresp()
        #         djSum.cmpbatch_lst = CmpBatchLst
        #         djSum.n_cmpbatches = NCmpBatches
        #         djSum.assay_id = AssayID
        #         NewEntry = True
        #     djSum.drval_type = 'MIC'
        #     djSum.act_types = row[('act_type','get_strList')]
        #     djSum.n_actives = row[('act_type','get_nAct')]
        #     djSum.act_score_ave = row[('act_score','mean')]
        #     djSum.inhibit_max_ave = row[('inhibit_max','mean')]
        #     djSum.pscore_ave = row[('pscore','mean')]
        #     djSum.drval_max    = row[('dr','DR_Range')]['Max']
        #     djSum.drval_min    = row[('dr','DR_Range')]['Min']
        #     djSum.drval_median = row[('dr','DR_Range')]['Median']
        #     djSum.n_assays = row[('dr','DR_Range')]['nDR']
        #     djSum.drval_unit   = row[('dr_unit','get_strList_unique')]

        #     if djSum.n_actives > 0:
        #         djCmpd.dr_n_actives += 1

        #     djSum.clean_Fields()
        #     validDict = djSum.validate()
        #     if validDict:
        #         validStatus = False
        #         # for k in validDict:
        #         #     print('Warning',k,validDict[k],'-')
        #         row.update(validDict)
        #         OutDict.append(row)

        #     if validStatus:
        #         if upload:
        #             if NewEntry or overwrite:
        #                 #djSum.chk_migration = 0
        #                 OutNumbers['Upload Entries'] += 1
        #                 djSum.save(user=appuser)    
    return(OutNumbers,OutDict)

# Summary DR from  AssayData MIC  ========================================================
def sum_cmpbatch_doseresponse_CC50(CmpBatchLst,NCmpBatches,djCmpd,
                                  upload=False,overwrite=False,appuser='J.Zuegg'):
# =========================================================================================
    validStatus = True
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)

    qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__result_type = 'MIC',
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values(
                                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
                                        'cc50','cc50_unit','act_type','act_score','pscore',
                                        'inhibit_max'
                                            )
    if qryCC50.exists():
    
        dfDR = pd.DataFrame(qryCC50).assign(cmpbatchs=CmpBatchs)
        dfDR.columns = ['plate_id','well_id','result_type','assay_id',
                        'dr','dr_unit','act_type','act_score','pscore',
                        'inhibit_max',
                        'cmpbatchs']

        OutDict = _process_sum_cmpbatch_dr('CC50',qryCC50,CmpBatchLst,NCmpBatches,CmpBatchs,djCmpd,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )

    return(OutNumbers,OutDict)

# Summary DR from  AssayData MIC  ========================================================
def sum_cmpbatch_doseresponse_HC50(CmpBatchLst,NCmpBatches,djCmpd,
                                  upload=False,overwrite=False,appuser='J.Zuegg'):
# =========================================================================================
    validStatus = True
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)

    qryCC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                    cmpbatch_lst__contains = CmpBatchLst, 
                                    n_cmpbatches = NCmpBatches, 
                                    testplate_id__result_type = 'MIC',
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values(
                                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
                                        'hc50','hc50_unit','act_type','act_score','pscore',
                                        'inhibit_max'
                                            )
    if qryCC50.exists():
    
        dfDR = pd.DataFrame(qryCC50).assign(cmpbatchs=CmpBatchs)
        dfDR.columns = ['plate_id','well_id','result_type','assay_id',
                        'dr','dr_unit','act_type','act_score','pscore',
                        'inhibit_max',
                        'cmpbatchs']

        OutDict = _process_sum_cmpbatch_dr('HC50',qryCC50,CmpBatchLst,NCmpBatches,CmpBatchs,djCmpd,OutNumbers,
                                            upload=upload,overwrite=overwrite,appuser=appuser )

    return(OutNumbers,OutDict)

# Summary DR Function  =======================================================================
def sum_cmpbatch_doseresponse(CmpBatchLst,upload=False,overwrite=False,appuser='J.Zuegg'):
# =========================================================================================
    OutNumbers = {'Processed':0,'New Entry':0, 'Upload Entries':0}
    OutDict = []

    CmpBatchs = COMPOUND_SEP.join(CmpBatchLst)
    NCmpBatches = len(CmpBatchLst)

    # Sum_Cmpd ---------------------------------------------------------------
    djCmpd = Summary_CmpBatch.get(CmpBatchLst,verbose=0)
    if djCmpd is None:
        djCmpd = Summary_CmpBatch()
        djCmpd.cmpbatch_lst = CmpBatchLst
        djCmpd.n_cmpbatches = len(CmpBatchLst)

    djCmpd.dr_n_assayids = 0
    djCmpd.dr_n_actives = 0

    # MIC
    _numbers,_outdict  = sum_cmpbatch_doseresponse_MIC(CmpBatchLst,NCmpBatches,djCmpd,
                                                upload=upload,overwrite=overwrite,appuser=appuser)
    if _outdict:
        OutDict = OutDict + _outdict
    for k in OutNumbers.keys():
        OutNumbers[k] += _numbers[k]

    # CC50
    _numbers,_outdict  = sum_cmpbatch_doseresponse_CC50(CmpBatchLst,NCmpBatches,djCmpd,
                                                upload=upload,overwrite=overwrite,appuser=appuser)
    if _outdict:
        OutDict = OutDict + _outdict
    for k in OutNumbers.keys():
        OutNumbers[k] += _numbers[k]

    # HC50
    _numbers,_outdict  = sum_cmpbatch_doseresponse_HC50(CmpBatchLst,NCmpBatches,djCmpd,
                                                upload=upload,overwrite=overwrite,appuser=appuser)
    if _outdict:
        OutDict = OutDict + _outdict
    for k in OutNumbers.keys():
        OutNumbers[k] += _numbers[k]

    # AssayData MIC ---------------------------------------------------------------
    # qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
    #                                 cmpbatch_lst__contains = CmpBatchLst, 
    #                                 n_cmpbatches = NCmpBatches, 
    #                                 testplate_id__result_type = 'MIC',
    #                                 testplate_id__plate_quality = 'Valid'                                            
    #                                 ).values(
    #                                     'testplate_id','testwell_id','testplate_id__result_type','testplate_id__assay_id',
    #                                     'mic','mic_unit','act_type','act_score','pscore',
    #                                     'inhibit_max'
    #                                         )
    # if qryMIC.exists():
    #     dfDR = pd.DataFrame(qryMIC).assign(cmpbatchs=CmpBatchs)
    #     dfDR.columns = ['plate_id','well_id','result_type','assay_id',
    #                     'mic','mic_unit','act_type','act_score','pscore',
    #                     'inhibit_max',
    #                     'cmpbatchs']

    #     pivDF = dfDR.groupby(['assay_id']).agg({'mic': [DR_Range],
    #                                         'inhibit_max': ['mean'],
    #                                         'pscore': ['mean'],
    #                                         'act_score': ['mean'],
    #                                         'act_type': [get_strList, get_nAct ],
    #                                         'mic_unit': [get_strList, get_strList_unique ],
    #                                         })
    #     #print( pivDF.columns)
    #     for AssayID,row in pivDF.iterrows():
    #         OutNumbers['Processed'] += 1
    #         djCmpd.dr_n_assayids += 1

    #         NewEntry = False
    #         djSum = Summary_CmpBatch_Doseresp.get(CmpBatchLst,AssayID,Exact=True,verbose=0)
    #         if djSum is None:
    #             djSum = Summary_CmpBatch_Doseresp()
    #             djSum.cmpbatch_lst = CmpBatchLst
    #             djSum.n_cmpbatches = NCmpBatches
    #             djSum.assay_id = AssayID
    #             NewEntry = True
    #         djSum.act_types = row[('act_type','get_strList')]
    #         djSum.n_actives = row[('act_type','get_nAct')]
    #         djSum.act_score_ave = row[('act_score','mean')]
    #         djSum.inhibit_max_ave = row[('inhibit_max','mean')]
    #         djSum.pscore_ave = row[('pscore','mean')]
    #         djSum.drval_type = 'MIC'
    #         djSum.drval_max    = row[('mic','DR_Range')]['Max']
    #         djSum.drval_min    = row[('mic','DR_Range')]['Min']
    #         djSum.drval_median = row[('mic','DR_Range')]['Median']
    #         djSum.n_assays = row[('mic','DR_Range')]['nDR']
    #         djSum.drval_unit   = row[('mic_unit','get_strList_unique')]

    #         # djSum, NewEntry = assign_sum_cmpbatch_doseresponse_MIC(row,CmpBatchLst,NCmpBatches,idx)

    #         if djSum.n_actives > 0:
    #             djCmpd.dr_n_actives += 1

    #         djSum.clean_Fields()
    #         validDict = djSum.validate()
    #         if validDict:
    #             validStatus = False
    #             # for k in validDict:
    #             #     print('Warning',k,validDict[k],'-')
    #             row.update(validDict)
    #             OutDict.append(row)

    #         if validStatus:
    #             if upload:
    #                 if NewEntry or overwrite:
    #                     #djSum.chk_migration = 0
    #                     OutNumbers['Upload Entries'] += 1
    #                     djSum.save(user=appuser)

    # Sum_Cmpd ---------------------------------------------------------------
    # djCmpd.dr_assayid_lst  = 
    # djCmpd.dr_actives_lst  = 
    if upload:
        djCmpd.save()

    return(OutNumbers,OutDict)