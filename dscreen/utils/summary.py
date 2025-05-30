#from dsample.models import CmpBatchList_Base, Project, COADD_Compound, ABase_Compound_Batch
from dplate.models import MasterPlate, TestPlate, TestWell
#from dchem.models import Chem_Structure
from dscreen.models import Screen_Run, AssayData_MIC, AssayData_CC50, AssayData_HC50
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
    if djRun.n_testplates > 0:
        djRun.screen_date = TestPlate.objects.filter(run_id = djRun.run_id).values('test_date').latest('test_date')['test_date']

