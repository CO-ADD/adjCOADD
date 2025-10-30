#from dsample.models import CmpBatchList_Base, Project, COADD_Compound, ABase_Compound_Batch
from dplate.models import MasterPlate, TestPlate, TestWell
#from dchem.models import Chem_Structure
from dscreen.models import Screen_Run, AssayData_MIC, AssayData_CC50, AssayData_HC50
from dsample.models import Project, Compound_Batch, COADD_Compound, ABase_Compound_Batch
from dgene.models import Genome_Sequence

#from applib.bio.bio_data import pScore, ActScore_DR, ActScore_SC

#from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

#--------------------------------------------------------------------
def update_screenrun_summary(djRun):
    """
     Update Calculated fields in Screen_Run
    """

    djRun.n_compounds = TestWell.objects.filter(plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
                                            ).values('cmpbatch_lst').distinct().count()
    # djRun.n_qc = 
    djRun.n_seq = Genome_Sequence.objects.filter(run_id = djRun.run_id).distinct().count()
    # djRun.n_structure = 
    #djRun.n_projects = 
    djRun.n_motherplates = MasterPlate.objects.filter(run_id=djRun.run_id).count()
    djRun.n_testplates = TestPlate.objects.filter(run_id=djRun.run_id).count()
    djRun.n_assays = TestPlate.objects.filter(run_id = djRun.run_id
                                            ).values('assay_id').distinct().count()
    djRun.n_inhibitions = TestWell.objects.filter(plate_id__result_type='Inhibition', 
                                            plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
                                            ).values('cmpbatch_lst').distinct().count()
    djRun.n_mic  = AssayData_MIC.objects.filter(run_id = djRun.run_id).count()
    djRun.n_cc50 = AssayData_CC50.objects.filter(run_id = djRun.run_id).count()
    djRun.n_hc50 = AssayData_HC50.objects.filter(run_id = djRun.run_id).count()
    
    # djRun.n_mic = TestWell.objects.filter(plate_id__result_type='MIC', 
    #                                         plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djRun.n_cc50 = TestWell.objects.filter(plate_id__result_type='CC50', 
    #                                         plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djRun.n_hc50 = TestWell.objects.filter(plate_id__result_type='HC50', 
    #                                         plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    
    djRun.n_synmic = TestWell.objects.filter(plate_id__result_type='synMIC', 
                                            plate_id__run_id = djRun.run_id, n_cmpbatches__gt = 0
                                            ).values('cmpbatch_lst').distinct().count()
    # Process Status of Testplates (Run)
    if str(djRun.run_type) in ['HCR','PSR']:
        # 0 No data
        # 1 Reads (Motherplates) -> Testplatelist
        # 2 Assays, Layout, MotherPlates -> Assign Cmpounds
        # 3 Compounds -> Calculate Inhibition
        # 4 Inhibition (Doseresponse)
        #  
        if djRun.n_testplates > 0:
            djRun.process_status = 1
        if djRun.n_assays > 0: 
            djRun.process_status = 2
        if djRun.n_compounds > 0: 
            djRun.process_status = 3
        if djRun.n_inhibitions > 0: 
            djRun.process_status = 4
            
    elif str(djRun.run_type) in ['SEQ']:
        if djRun.n_seq > 0:
            djRun.process_status = 1
    
    if djRun.n_testplates > 0:
        djRun.screen_date = TestPlate.objects.filter(run_id = djRun.run_id).values('test_date').latest('test_date')['test_date']


#--------------------------------------------------------------------
def get_projects_screenrun(djRun):
#--------------------------------------------------------------------
    #
    # Return a list of dictionary of Projects in ScreenRun
    #   [{'project_id','project_type'}]
    #    
    # from COADD_Compound
    # 
    _projects = []
    _cmpbatch_lsts = TestWell.objects.filter(plate_id__run_id = djRun.run_id, 
                                             n_cmpbatches__gt = 0).values('cmpbatch_lst').distinct().order_by()
    
    _cmp_lists = {'COADD':[],'ABASE':[],'LIBRARY':[]}
    for _cb in _cmpbatch_lsts:
        _cmpbatch_lst = _cb['cmpbatch_lst']
        for _c in _cb['cmpbatch_lst']:
            djBatch = Compound_Batch.get(_c)
            _cmp_lists[djBatch.batch_source].append(_c)
    
    if len(_cmp_lists['COADD'])>0:
        _prj_lst = COADD_Compound.objects.filter(compound_id__in=_cmp_lists['COADD']).values('project_id').distinct().order_by()
        for _p in _prj_lst:
            _projects.append({'project_id':_p['project_id'],'project_type':'CO-ADD'})
             
    # if len(_cmp_lists['ABASE'])>0:
    #     _prj_lst = ABase_Compound_Batch.objects.filter(cmpbatch_id=_cmp_lists['ABASE']).values('project_id').distinct().order_by()
    #     for _p in _prj_lst:
    #         _projects.append({'project_id':_p['project_id'],'project_type':'ABASE'})
            
    return(_projects)