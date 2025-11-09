from dplate.models import MasterPlate, TestPlate, TestWell, MasterWell
from dsample.models import Compound_Batch, COADD_Compound, ABase_Compound_Batch

#------------------------------------------------
def update_project_summary(djPrj):
    """
     Update Calculated fields in Project
    """
    
    _COADD_Compounds = COADD_Compound.objects.filter(project_id = djPrj.project_id)
    #COADD_Compound.objects.filter(project_id = djPrj.project_id).count()
    
    djPrj.n_compounds = _COADD_Compounds.count()
    djPrj.n_barcodes = 0
    
    if djPrj.n_compounds > 0 :
        _CmpBatches = Compound_Batch.objects.filter(cmpbatch_id__in = _COADD_Compounds)
        _Barcodes = MasterWell.objects.filter(cmpbatch_id__in = _CmpBatches, barcode__isnull=False).values_list('barcode')
        djPrj.n_barcodes += _Barcodes.count()

    # djPrj.n_mcc_compounds = ABase_Compound_Batch.objects.filter(project_id = djPrj.project_id).count()
    #     _CmpBatches = Compound_Batch.objects.filter(cmpbatch_id__in = ABase_Compound_Batch)
    #     _Barcodes = MasterWell.objects.filter(cmpbatch_id__in = _CmpBatches, plate_id__plate_type = 'Storage')
    #     djPrj.n_barcodes += _Barcodes.count()

    

    # djPrj.n_qc = 
    # djPrj.n_structure = 
    #djPrj.n_projects = 
    # djPrj.n_motherplates = MasterPlate.objects.filter(run_id=djPrj.run_id).count()
    # djPrj.n_testplates = TestPlate.objects.filter(run_id=djPrj.run_id).count()
    # djPrj.n_assays = TestPlate.objects.filter(run_id = djPrj.run_id
    #                                         ).values('assay_id').distinct().count()
    # djPrj.n_inhibitions = TestWell.objects.filter(plate_id__result_type='Inhibition', 
    #                                         plate_id__run_id = djPrj.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djPrj.n_mic  = AssayData_MIC.objects.filter(run_id = djPrj.run_id).count()
    # djPrj.n_cc50 = AssayData_CC50.objects.filter(run_id = djPrj.run_id).count()
    # djPrj.n_hc50 = AssayData_HC50.objects.filter(run_id = djPrj.run_id).count()
    
    # djPrj.n_mic = TestWell.objects.filter(plate_id__result_type='MIC', 
    #                                         plate_id__run_id = djPrj.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djPrj.n_cc50 = TestWell.objects.filter(plate_id__result_type='CC50', 
    #                                         plate_id__run_id = djPrj.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djPrj.n_hc50 = TestWell.objects.filter(plate_id__result_type='HC50', 
    #                                         plate_id__run_id = djPrj.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    
    # djPrj.n_synmic = TestWell.objects.filter(plate_id__result_type='synMIC', 
    #                                         plate_id__run_id = djPrj.run_id, n_cmpbatches__gt = 0
    #                                         ).values('cmpbatch_lst').distinct().count()
    # djPrj.screen_date = TestPlate.objects.filter(run_id = djPrj.run_id).values('test_date').latest('test_date')['test_date']

