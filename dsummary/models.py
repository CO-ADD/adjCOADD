from django.db import models
import re

#from django_rdkit import models
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from applib.data.str_lists import strList_to_List
#from dchem.models import Chem_Structure
from dsample.models import CmpBatchList_Base, Project, COADD_Compound, ABase_Compound_Batch
from dplate.models import MasterPlate, TestPlate, TestWell
from dchem.models import Chem_Structure
from dscreen.models import Screen_Run, AssayData_MIC, AssayData_CC50, AssayData_HC50
from applib.bio.bio_data import pScore, ActScore_DR, ActScore_SC

from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

      
#-------------------------------------------------------------------------------------------------
# Summary Project
#-------------------------------------------------------------------------------------------------

# class Summary_Project(AuditModel):
#     """
#     List of Summary for each ScreenRun
#     """

#     #---------------------------------------------------------------------------------------------
#     project_id = models.OneToOneField(Project, primary_key=True, verbose_name = "Project ID", on_delete=models.DO_NOTHING,
#                                         db_column="project_id", related_name="%(class)s_project_id")

#     n_compounds = models.IntegerField(default=0, verbose_name = "#Cpmds")
#     n_mcc_compounds = models.IntegerField(default=0, verbose_name = "#MCC")
#     n_barcode = models.IntegerField(default=0, verbose_name = "#BCode")
#     n_structure = models.IntegerField(default=0, verbose_name = "#Struc")
#     n_motherplates = models.IntegerField(default=0, verbose_name = "#MP")
#     # n_testplates = models.IntegerField(default=0, verbose_name = "#TP")
#     n_runids = models.IntegerField(default=0, verbose_name = "#Runs")
#     n_assays = models.IntegerField(default=0, verbose_name = "#Assays")
#     n_ps_compounds = models.IntegerField(default=0, verbose_name = "#PS")
#     n_dr_compounds = models.IntegerField(default=0, verbose_name = "#DR")
#     n_syn_compounds = models.IntegerField(default=0, verbose_name = "#SYN")
#     screen_date = models.DateField(null=True, blank=True, verbose_name="Screen Date")
#     n_sc_hits = models.IntegerField(default=0, verbose_name = "#Inhib Hits")
#     n_mic_hits = models.IntegerField(default=0, verbose_name = "#MIC Hits")
#     n_tox_hits = models.IntegerField(default=0, verbose_name = "#Tox Hits")

#     #------------------------------------------------
#     class Meta:
#         app_label = 'dsummary'
#         db_table = 'sum_project'
#         ordering=['project_id']
#         indexes = [
#             models.Index(name="sprj_ncmp_idx", fields=['n_compounds']),
#             models.Index(name="sprj_nstr_idx", fields=['n_structure']),
#             models.Index(name="sprj_nsh_idx", fields=['n_sc_hits']),
#             models.Index(name="sprj_nmh_idx", fields=['n_mic_hits']),
#             models.Index(name="sprj_nth_idx", fields=['n_tox_hits']),
#             # models.Index(name="scmpsc_ascr_idx", fields=['act_score_ave']),
#             # models.Index(name="scmpsc_inhin_idx", fields=['inhibition_ave']),
#             # models.Index(name="scmpsc_mscr_idx", fields=['mscore_ave']),
#         ]

#     #------------------------------------------------
#     def __str__(self) -> str:
#         return f"{self.project_id}"

#     #------------------------------------------------
#     def __repr__(self) -> str:
#         return f"{self.project_id} {self.n_compounds}"

#     #------------------------------------------------
#     @classmethod
#     def update(cls,ProjectID, to_save=False,verbose=0):
#         retObj = cls.get(ProjectID)
#         if retObj is None:
#             retObj = cls()
#             retObj.project_id = Project.get(ProjectID)
            
#         if retObj.project_id:
#             retObj.update_summary()            
#             if to_save:
#                 retObj.save()
        
#     #------------------------------------------------
#     def update_summary(self):

#         self.n_compounds = COADD_Compound.objects.filter(project_id = self.project_id).count()
#         self.n_mcc_compounds = ABase_Compound_Batch.objects.filter(project_id = self.project_id).count()
#         # self.n_qc = 
#         # self.n_structure = 
#         #self.n_projects = 
#         # self.n_motherplates = MasterPlate.objects.filter(run_id=self.run_id).count()
#         # self.n_testplates = TestPlate.objects.filter(run_id=self.run_id).count()
#         # self.n_assays = TestPlate.objects.filter(run_id = self.run_id
#         #                                         ).values('assay_id').distinct().count()
#         # self.n_inhibitions = TestWell.objects.filter(plate_id__result_type='Inhibition', 
#         #                                         plate_id__run_id = self.run_id, n_cmpbatches__gt = 0
#         #                                         ).values('cmpbatch_lst').distinct().count()
#         # self.n_mic  = AssayData_MIC.objects.filter(run_id = self.run_id).count()
#         # self.n_cc50 = AssayData_CC50.objects.filter(run_id = self.run_id).count()
#         # self.n_hc50 = AssayData_HC50.objects.filter(run_id = self.run_id).count()
        
#         # self.n_mic = TestWell.objects.filter(plate_id__result_type='MIC', 
#         #                                         plate_id__run_id = self.run_id, n_cmpbatches__gt = 0
#         #                                         ).values('cmpbatch_lst').distinct().count()
#         # self.n_cc50 = TestWell.objects.filter(plate_id__result_type='CC50', 
#         #                                         plate_id__run_id = self.run_id, n_cmpbatches__gt = 0
#         #                                         ).values('cmpbatch_lst').distinct().count()
#         # self.n_hc50 = TestWell.objects.filter(plate_id__result_type='HC50', 
#         #                                         plate_id__run_id = self.run_id, n_cmpbatches__gt = 0
#         #                                         ).values('cmpbatch_lst').distinct().count()
        
#         # self.n_synmic = TestWell.objects.filter(plate_id__result_type='synMIC', 
#         #                                         plate_id__run_id = self.run_id, n_cmpbatches__gt = 0
#         #                                         ).values('cmpbatch_lst').distinct().count()
#         # self.screen_date = TestPlate.objects.filter(run_id = self.run_id).values('test_date').latest('test_date')['test_date']

#-------------------------------------------------------------------------------------------------
# Summary Screening Data Models
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
class Summary_CmpBatch(CmpBatchList_Base):
    """
    List of Summary Activity for each CmpBatch
    """
    ASSAY_CLASSES = {'gp':0,'gn':1,'fg':2,'cc':3,'hc':4,'gnm':5}

    GNM_ASSAYS = ['GN_0046','GN_0048','GN_0049','GN_0211']

    sc_n_assayids = models.SmallIntegerField(default=-1, blank=True, verbose_name = "SC #AssayIDs")
    sc_n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "SC #Actives")

    sc_assayid_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "SC AssayID List")
    sc_actives_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "SC Actives List")

    dr_n_assayids = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#AssayIDs")
    dr_n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    dr_assayid_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "DR AssayID List")
    dr_actives_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "DR Actives List")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_cmpbatch'
        indexes = [
            GinIndex(name="scmp_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="scmp_scassid_idx", fields=['sc_n_assayids']),
            models.Index(name="scmp_scnact_idx", fields=['sc_n_actives']),
            models.Index(name="scmp_drassid_idx", fields=['dr_n_assayids']),
            models.Index(name="scmp_drnact_idx", fields=['dr_n_actives']),
            # models.Index(name="scmpsc_ascr_idx", fields=['act_score_ave']),
            # models.Index(name="scmpsc_inhin_idx", fields=['inhibition_ave']),
            # models.Index(name="scmpsc_mscr_idx", fields=['mscore_ave']),
        ]


#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
class Summary_CmpBatch_Inhib(CmpBatchList_Base):
    """
    List of Summary Activity for each CmpBatch
    """
#-------------------------------------------------------------------------------------------------
 
    # Assay Conditions
    sum_assay_id = models.CharField(max_length=100, blank=True, verbose_name = "SumAssay ID")
    n_assays = models.IntegerField(default=-1, blank=True, verbose_name = "#Assay")

    # Activity Summary
    act_types = models.CharField(max_length=250, blank=True, verbose_name = "Active Tupes")
    # active_lst
    n_actives = models.IntegerField(default=-1, blank=True, verbose_name = "#Actives")
    act_score = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "Act Score")

    inhibition_ave = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Ave")
    inhibition_std = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Std")
    inhibition_min = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Min")
    inhibition_max = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Max")
    # inhibition_lst
    # inhibition_range = models.CharField(max_length=50, blank=True, verbose_name = "Inhibition Range")
    mscore_ave = models.DecimalField(max_digits=9, decimal_places=3, verbose_name = "MScore Max")

    concs_lst = models.CharField(max_length=250, blank=True, verbose_name = "Concs")
    # conc_unit_lst

    # Summary Meta data
    # run_id_lst ArrayField(models.CharField(max_length=15, default="", db_index = True), 
    #                             size=MAX_CMPBATCHES, verbose_name = "CmpBatch List", null=True, blank=True)
    # run_id_date = models.DateField(null=True, blank=True, verbose_name = "Date")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_cmpbatch_sc'
        ordering=['sum_assay_id']
        indexes = [
            GinIndex(name="scmpsc_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="scmpsc_assid_idx", fields=['sum_assay_id']),
            models.Index(name="scmpsc_nact_idx", fields=['n_actives']),
            models.Index(name="scmpsc_ascr_idx", fields=['act_score']),
            models.Index(name="scmpsc_inhin_idx", fields=['inhibition_ave']),
            models.Index(name="scmpsc_mscr_idx", fields=['mscore_ave']),
        ]

    #------------------------------------------------
    # Returns an User instance if found by name
    #------------------------------------------------
    @classmethod
    def get(cls,CmpBatchLst,AssayID, Exact=True, verbose=0):
        try:
            if Exact:
                retInstance = cls.objects.get(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst), sum_assay_id=AssayID)
            else:
                retInstance = cls.objects.get(cmpbatch_lst__contains=CmpBatchLst, sum_assay_id=AssayID)
        except:
            if verbose:
                logger.warning(f"[Summary CmpBatch Inhibition Not Found] {CmpBatchLst} {AssayID}")
            retInstance = None
        return(retInstance)


    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,CmpBatchLst,AssayID, Exact=True, verbose=0):
        if Exact:
            return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst), assay_id=AssayID).exists()
        else:
            return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, assay_id=AssayID).exists()

    #------------------------------------------------
    # Sets the Act_Score
    def set_actscores(self,verbose=0):
        self.act_score = ActScore_SC(self.inhibition_ave,ZScore=self.mscore_ave)

#-------------------------------------------------------------------------------------------------
class Summary_CmpBatch_Doseresp(CmpBatchList_Base):
    """
    List of Summary Activity for each CmpBatch
    """
#-------------------------------------------------------------------------------------------------

    # Assay Conditions
    sum_assay_id = models.CharField(max_length=100, blank=True, verbose_name = "SumAssay ID")
    n_assays = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Assay")

    # Activity Summary
    act_types = models.CharField(max_length=250, blank=True, verbose_name = "Active Types")
    # active_lst
    n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    act_score = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "Act Score Ave")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore Ave")

    inhibit_max_ave = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Max Ave")
    #inhibit_maxs  = models.CharField(max_length=1024, blank=False, verbose_name = "DRs")
    drval_type    = models.CharField(max_length=15, blank=False, verbose_name = "DR Type")
    drval_max    = models.CharField(max_length=20, blank=False, verbose_name = "DR Max")
    drval_min    = models.CharField(max_length=20, blank=False, verbose_name = "DR Min")
    drval_median = models.CharField(max_length=20, blank=False, verbose_name = "DR Median")
    drval_unit   = models.CharField(max_length=25, blank=False, verbose_name = "DR Unit")

    drval_std_geomean = models.CharField(max_length=20, blank=False, verbose_name = "DR Std Geomean")
    drval_std_unit    = models.CharField(max_length=25, blank=False, verbose_name = "DR Std Unit")

    # Summary Meta data
    # run_id_lst ArrayField(models.CharField(max_length=15, default="", db_index = True), 
    #                             size=MAX_CMPBATCHES, verbose_name = "CmpBatch List", null=True, blank=True)
    # run_id_date = models.DateField(null=True, blank=True, verbose_name = "Date")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_cmpbatch_dr'
        ordering=['sum_assay_id']
        indexes = [
            GinIndex(name="scmpdr_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="scmpdr_assid_idx", fields=['sum_assay_id']),
            models.Index(name="scmpdr_nact_idx", fields=['n_actives']),
            models.Index(name="scmpdr_ascr_idx", fields=['act_score']),
            models.Index(name="scmpdr_pscr_idx", fields=['pscore']),
            models.Index(name="scmpdr_drt_idx", fields=['drval_type']),
            # models.Index(name="scmpdr_mscr_idx", fields=['mscore_ave']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.id}"
        

    #------------------------------------------------
    # Returns an User instance if found by name
    #------------------------------------------------
    @classmethod
    def get(cls,CmpBatchLst, AssayID, Exact=True, verbose=0):
        try:
            if Exact:
                retInstance = cls.objects.get(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst), sum_assay_id=AssayID)
            else:
                retInstance = cls.objects.get(cmpbatch_lst__contains=CmpBatchLst, sum_assay_id=AssayID)
        except:
            if verbose:
                logger.warning(f"[Summary CmpBatch DoseResp Not Found] {CmpBatchLst} {AssayID}")
            retInstance = None
        return(retInstance)


    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,CmpBatchLst,AssayID, Exact=True, verbose=0):
        if Exact:
            return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst), assay_id=AssayID).exists()
        else:
            return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, assay_id=AssayID).exists()

    #------------------------------------------------
    # Sets the Act_Score and pScores
    def set_actscores(self,verbose=0):
        if self.drval_type == 'HC50':
            self.act_score = ActScore_DR(self.drval_median,self.drval_unit,DMax=self.inhibit_max_ave,cutoff_inhib=10)
        else:
            self.act_score = ActScore_DR(self.drval_median,self.drval_unit,DMax=self.inhibit_max_ave)
        self.pscore = pScore(self.drval_std_geomean,self.drval_std_unit,self.inhibit_max_ave,MW=self.cmpbatch_id.full_mw,gtShift=3,drMax2=40)

#-------------------------------------------------------------------------------------------------
class Summary_Structure(AuditModel):
    """
    List of Summary Activity for each CmpBatch
    """

    ASSAY_CLASSES = ['GP','GN','GNMemb','FG','CC','HC']

    #---------------------------------------------------------------------------------------------
    structure_id = models.OneToOneField(Chem_Structure, primary_key=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
                                        db_column="structure_id", related_name="%(class)s_structure_id")
    sc_n_assayids = models.SmallIntegerField(default=-1, blank=True, verbose_name = "SC #AssayIDs")
    sc_n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "SC #Actives")

    sc_assayid_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "SC AssayID List")
    sc_actives_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "SC Actives List")

    dr_n_assayids = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#AssayIDs")
    dr_n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    dr_assayid_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "DR AssayID List")
    dr_actives_lst  = ArrayField(models.SmallIntegerField(default=-1, blank=True),
                                size=len(ASSAY_CLASSES), null=True, blank=True, db_index=True, verbose_name = "DR Actives List")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_structure'
        ordering=['structure_id']
        indexes = [
            models.Index(name="sstr_scassid_idx", fields=['sc_n_assayids']),
            models.Index(name="sstr_scnact_idx", fields=['sc_n_actives']),
            models.Index(name="sstr_drassid_idx", fields=['dr_n_assayids']),
            models.Index(name="sstr_drnact_idx", fields=['dr_n_actives']),
            # models.Index(name="scmpsc_ascr_idx", fields=['act_score_ave']),
            # models.Index(name="scmpsc_inhin_idx", fields=['inhibition_ave']),
            # models.Index(name="scmpsc_mscr_idx", fields=['mscore_ave']),
        ]

    #------------------------------------------------
    # Returns an User instance if found by name
    #------------------------------------------------
    @classmethod
    def get(cls,StructureID, verbose=0):
        try:
            retInstance = cls.objects.get(structure_id=StructureID)
        except:
            if verbose:
                logger.warning(f"[Summary Structure Inhibition Not Found] {StructureID} ")
            retInstance = None
        return(retInstance)

#-------------------------------------------------------------------------------------------------
class Summary_Structure_Inhib(AuditModel):
    """
    List of Summary Activity for each Structure
    """
#-------------------------------------------------------------------------------------------------
    #from dplate.models import TestPlate

    # Structure ID
    structure_id= models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")

    # Assay Conditions
    sum_assay_id = models.CharField(max_length=25, blank=True, verbose_name = "SumAssay ID")
    n_assays = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Assay")

    # Activity Summary
    act_types = models.CharField(max_length=250, blank=True, verbose_name = "Active Tupes")
    # active_lst
    n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    act_score = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "Act Score")

    inhibition_ave = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Ave")
    inhibition_std = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Std")
    inhibition_min = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Min")
    inhibition_max = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Max")
    # inhibition_lst
    # inhibition_range = models.CharField(max_length=50, blank=True, verbose_name = "Inhibition Range")
    mscore_ave = models.DecimalField(max_digits=9, decimal_places=3, verbose_name = "MScore Max")

    # conc_lst
    # conc_unit_lst

    # Summary Meta data
    # run_id_lst ArrayField(models.CharField(max_length=15, default="", db_index = True), 
    #                             size=MAX_CMPBATCHES, verbose_name = "CmpBatch List", null=True, blank=True)
    # run_id_date = models.DateField(null=True, blank=True, verbose_name = "Date")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_structure_sc'
        ordering=['sum_assay_id','structure_id']
        constraints = [
            models.UniqueConstraint(name='sstrsc_cst', fields=['structure_id', 'sum_assay_id'], )
        ]        
        indexes = [
#            GinIndex(name="scmpsc_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="sstrsc_sid_idx", fields=['structure_id']),
            models.Index(name="sstrsc_assid_idx", fields=['sum_assay_id']),
            models.Index(name="sstrsc_nact_idx", fields=['n_actives']),
            models.Index(name="sstrsc_ascr_idx", fields=['act_score']),
            models.Index(name="sstrsc_inhin_idx", fields=['inhibition_ave']),
            models.Index(name="sstrsc_mscr_idx", fields=['mscore_ave']),
        ]

    #------------------------------------------------
    # Returns an User instance if found by name
    #------------------------------------------------
    @classmethod
    def get(cls,StructureID,AssayID, verbose=0):
        try:
            retInstance = cls.objects.get(structure_id=StructureID, sum_assay_id=AssayID)
        except:
            if verbose:
                logger.warning(f"[Summary Structure Inhibition Not Found] {StructureID} {AssayID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,StructureID,AssayID, verbose=0):
        return cls.objects.filter(structure_id=StructureID, assay_id=AssayID).exists()

    #------------------------------------------------
    # Sets the Act_Score
    def set_actscores(self,verbose=0):
        self.act_score = ActScore_SC(self.inhibition_ave,ZScore=self.mscore_ave)


#-------------------------------------------------------------------------------------------------
class Summary_Structure_Doseresp(AuditModel):
    """
    List of Summary Activity for each Structure
    """
#-------------------------------------------------------------------------------------------------

    # Structure ID
    structure_id= models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")

    # Assay Conditions
    sum_assay_id = models.CharField(max_length=100, blank=True, verbose_name = "SumAssay ID")
    n_assays = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Assay")

    # Activity Summary
    act_types = models.CharField(max_length=250, blank=True, verbose_name = "Active Types")
    # active_lst
    n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    act_score = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    inhibit_max_ave = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Max Ave")
    #inhibit_maxs  = models.CharField(max_length=1024, blank=False, verbose_name = "DRs")
    drval_type    = models.CharField(max_length=15, blank=False, verbose_name = "DR Type")

    drval_max    = models.CharField(max_length=20, blank=False, verbose_name = "DR Max")
    drval_min    = models.CharField(max_length=20, blank=False, verbose_name = "DR Min")
    drval_median = models.CharField(max_length=20, blank=False, verbose_name = "DR Median")
    drval_unit   = models.CharField(max_length=25, blank=False, verbose_name = "DR Unit")

    drval_std_geomean = models.CharField(max_length=20, blank=False, verbose_name = "DR Std Geomean")
    drval_std_unit    = models.CharField(max_length=25, blank=False, verbose_name = "DR Std Unit")

    # Summary Meta data
    # run_id_lst ArrayField(models.CharField(max_length=15, default="", db_index = True), 
    #                             size=MAX_CMPBATCHES, verbose_name = "CmpBatch List", null=True, blank=True)
    # run_id_date = models.DateField(null=True, blank=True, verbose_name = "Date")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_structure_dr'
        ordering=['sum_assay_id','structure_id']
        constraints = [
            models.UniqueConstraint(name='sstrdr_cst', fields=['structure_id', 'sum_assay_id'], )
        ]        
        indexes = [
#            GinIndex(name="scmpsc_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="sstrdr_sid_idx", fields=['structure_id']),
            models.Index(name="sstrdr_assid_idx", fields=['sum_assay_id']),
            models.Index(name="sstrdr_nact_idx", fields=['n_actives']),
            models.Index(name="sstrdr_ascr_idx", fields=['act_score']),
            models.Index(name="sstrdr_pscr_idx", fields=['pscore']),
            models.Index(name="sstrdr_drt_idx", fields=['drval_type']),
        ]

    #------------------------------------------------
    # Returns an User instance if found by name
    #------------------------------------------------
    @classmethod
    def get(cls,StructureID,AssayID, verbose=0):
        try:
            retInstance = cls.objects.get(structure_id=StructureID, sum_assay_id=AssayID)
        except:
            if verbose:
                logger.warning(f"[Summary Structure DoseResponse Not Found] {StructureID} {AssayID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,StructureID,AssayID, verbose=0):
        return cls.objects.filter(structure_id=StructureID, assay_id=AssayID).exists()

    #------------------------------------------------
    # Sets the Act_Score and pScores
    def set_actscores(self,verbose=0):
        if self.drval_type == 'HC50':
            self.act_score = ActScore_DR(self.drval_median,self.drval_unit,DMax=self.inhibit_max_ave,cutoff_inhib=10)
        else:
            self.act_score = ActScore_DR(self.drval_median,self.drval_unit,DMax=self.inhibit_max_ave)
        self.pscore = pScore(self.drval_std_geomean,self.drval_std_unit,self.inhibit_max_ave,MW=self.structure_id.mw,gtShift=3,drMax2=40)
