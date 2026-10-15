
#from django_rdkit import models
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
from applib.data.str_lists import strList_to_List, split_StrList
from applib.bio.bio_data import pScore, ActScore_DR, ActScore_SC
from dsample.models import CmpBatchList_Base, Compound_Batch

from dcell.models import Cell
from dorganism.models import Organism
from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------------------------
# Screening Application Model
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
class Screen_Run(AuditModel):
    """
    List of Screening runs
    """
#-------------------------------------------------------------------------------------------------
    LIST_VIEW_FIELDS = {
        "run_id":{'Run ID': {'run_id':URL_LINKS['screenrun_id']}},
        "run_type":"Run Type",
        "assay_note":"Assay",
        "run_status":"Status",
        "run_project":"Project",
        "run_name":"Name",
        #"run_date":"Run Date",
        "run_conditions":"Conditions",
        "run_issues":"Issues",
        #"process_status":"Process Status",
        # Calculated Fields
        "screen_date": "Screen Date",
        "n_compounds":"#Cmpds",     
        "n_motherplates":"#MP",     
        "n_testplates":"#TP",     
        "n_testplates_valid": "#TP Valid",
        "n_assays":"#Ass",     
        "n_inhibitions":"#Inhib",     
        "n_mic":"#MIC",     
        "n_cc50":"#CC50",     
        "n_hc50":"#HC50",
        "n_synmic":"#micSyn",     
        "n_seq":"#Seq",     
    }

    DICTIONARY_FIELDS = {
        'run_type':'Run_Type',
        'run_status':'Process_Status',
    }

    VIEW_GROUPS = [
        ['run_type','run_status','run_name','run_source','run_project','run_date'],
        ['run_conditions','assay_note','run_issues'],
        ['n_compounds', 'n_structure','n_motherplates','n_testplates','n_testplates_valid','n_qc','n_seq'],
        ['screen_date','n_assays','n_inhibitions','n_mic','n_cc50','n_hc50','n_synmic','process_status']
    ]


    run_id = models.CharField(max_length=15, primary_key=True, blank=True, verbose_name = "Run ID")
    run_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Type", on_delete=models.DO_NOTHING,
        db_column="run_type", related_name="%(class)s_RunType+")
    run_name = models.CharField(max_length=500, blank=True, verbose_name = "Run Name")
    #run_folder = models.CharField(max_length=120, blank=True, verbose_name = "Run Name")
    assay_note = models.CharField(max_length=250, blank=True, verbose_name = "Assay Note")
    run_conditions = models.CharField(max_length=250, blank=True, verbose_name = "Run Conditions")
    run_issues = models.CharField(max_length=250, blank=True, verbose_name = "Run Issues")
    run_date = models.DateField(null=True, blank=True, verbose_name = "Run Date")
    run_source = models.CharField(max_length=50, blank=True, verbose_name = "Source")
    run_project = models.CharField(max_length=50, blank=True, verbose_name = "Project")
    run_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Status", on_delete=models.DO_NOTHING,
        db_column="run_status", related_name="%(class)s_run_status+")

    # Add PlatePrep FileName 

    # -- Calculated Fields - No View/Update ------------------------------------------
    CALCULATED_FIELDS = ['n_compounds','n_qc','n_structure',
                         'n_motherplates','n_testplates','n_testplates_valid','n_assays',
                         'n_inhibitions','n_mic','n_cc50','n_hc50','n_synmic',
                         'n_seq',
                         'screen_date',
                         'process_status']

    n_compounds = models.IntegerField(default=0, blank=True, verbose_name = "#Cpmds")
    #n_projects = models.SmallIntegerField(default=0, verbose_name = "#Projects")
    n_seq = models.IntegerField(default=0, blank=True, verbose_name = "#Seq")
    n_qc = models.IntegerField(default=0, blank=True, verbose_name = "#QC")
    n_structure = models.IntegerField(default=0, blank=True, verbose_name = "#Struc")
    n_motherplates = models.IntegerField(default=0, blank=True, verbose_name = "#MP")
    n_testplates = models.IntegerField(default=0, blank=True, verbose_name = "#TP")
    n_testplates_valid = models.IntegerField(default=0, blank=True, verbose_name = "#TP Valid")
    n_assays = models.IntegerField(default=0, blank=True, verbose_name = "#Assays")
    n_inhibitions = models.IntegerField(default=0, blank=True, verbose_name = "#Inhib")
    n_mic = models.IntegerField(default=0, blank=True, verbose_name = "#MIC")
    n_cc50 = models.IntegerField(default=0, blank=True, verbose_name = "#CC50")
    n_hc50 = models.IntegerField(default=0, blank=True, verbose_name = "#HC50")
    n_synmic = models.IntegerField(default=0, blank=True, verbose_name = "#micSyn")
    screen_date = models.DateField(null=True, blank=True, verbose_name="Screen Date")
    process_status = models.IntegerField(default=0, blank=True, verbose_name = "Process Status")


    #------------------------------------------------
    class Meta:
        app_label = 'dscreen'
        db_table = 'screen_run'
        ordering=['run_type','run_id']
        indexes = [
            models.Index(name="run_type_idx", fields=['run_type']),
            models.Index(name="run_ncp_idx", fields=['n_compounds']),
            models.Index(name="run_nst_idx", fields=['n_structure']),
            models.Index(name="run_nas_idx", fields=['n_assays']),
            models.Index(name="run_nmp_idx", fields=['n_motherplates']),
            models.Index(name="run_ntp_idx", fields=['n_testplates']),
            models.Index(name="run_nvtp_idx", fields=['n_testplates_valid']),
            models.Index(name="run_ninh_idx", fields=['n_inhibitions']),
            models.Index(name="run_pst_idx", fields=['process_status']),
        ]

    # #------------------------------------------------
    # def __str__(self) -> str:
    #     return f"{self.run_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.run_id} [{self.run_type}]"

    #------------------------------------------------
    @classmethod
    def get(cls,RunID,verbose=0):
        try:
            retInstance = cls.objects.get(run_id=RunID)
        except:
            if verbose:
                print(f"[Run_ID Not Found] {RunID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,RunID,verbose=0):
        return cls.objects.filter(run_id=RunID.strip()).exists()
    

    #------------------------------------------------
    @classmethod
    def str_RunID(cls,RunType,RunNo) -> str:
    #
    # Input:    RunClass PSR, HCR, QCR,...
    #           RunNo 
    # Output:   Run_ID as string like PSR00001 
    #
        return(f"{RunType}{RUN_SEP}{RunNo:05d}")

    #------------------------------------------------
    @classmethod
    def find_Next_RunID(cls,RunType,RunClassTypes = RUN_CLASSES) -> str:
        if RunType in RunClassTypes:
            Run_IDSq=Sequence(RunType)
            Run_nextID = next(Run_IDSq)
            Run_strID = cls.str_RunID(RunType,Run_nextID)
            while cls.exists(Run_strID):
                Run_nextID = next(Run_IDSq)
                Run_strID = cls.str_RunID(RunType,Run_nextID)
            return(Run_strID)    
        else:
            return(None)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        verbose = kwargs.get('verbose',0)

        if not self.run_id:
            self.run_id = self.find_Next_RunID(str(self.run_type.dict_value))
        if self.run_id:
            if verbose>0:
                print(f'Saving ScreenRun [{self.run_id}] [{self.run_type.dict_value}]')
            super(Screen_Run, self).save(*args, **kwargs)
        # else:
        #     super(Screen_Run, self).save(*args, **kwargs) 


#-------------------------------------------------------------------------------------------------
class Assay(AuditModel):
    """
    List of Assays
    """
#-------------------------------------------------------------------------------------------------
    LIST_VIEW_FIELDS = {
        "assay_id":"Assay ID",
        "assay_code":"Assay Code",
        #"assay_subtype":"Assay SubType",
        "organism_id":"Organism",
        "cell_id":"Cell",
        # "run_project":"Project",
        'test_media' :"Media",
        'test_dye' :"Dye/Kit",
        'test_enviroment' : "Enviroment",
        'test_temperature' : "Temp",
        'test_time' :  "Test Time",
        'test_additive' : "Additive",
        'subculture_type' : "Subculture",
        #'incubation_time' : "Incubation Time",
    }

    DICTIONARY_FIELDS = {
        # 'run_type':'Run_Type',
        # 'run_status':'Process_Status',
    }

    CALCULATED_FIELDS = []
    VIEW_GROUPS=[]

    assay_id = models.CharField(max_length=100,primary_key=True, verbose_name = "Assay ID")
    ora_assay_id = models.CharField(max_length=100,blank=True, verbose_name = "Ora Assay ID")
    assay_type = models.CharField(max_length=30, verbose_name = "AssayType" )
    assay_subtype = models.CharField(max_length=50, verbose_name = "AssaySubType" )
    assay_panel = ArrayField(models.CharField(max_length=20, blank=True),size=30, null=True, blank=True, verbose_name = "Assay Panel")
    sum_assay_id  = models.CharField(max_length=100, blank=True, verbose_name = "Assay ID for Summary")
    coadd_assay_id  = models.CharField(max_length=15, blank=True, verbose_name = "COADD Assay ID")
    assay_note = models.CharField(max_length=150, blank=True, verbose_name = "Assay Note")
    assay_code = models.CharField(max_length=20, blank=True, verbose_name = "Assay Code")
    test_media = models.CharField(max_length=150, blank=True, verbose_name = "Media")
    test_dye = models.CharField(max_length=150, blank=True, verbose_name = "Dye/Kit")
    test_enviroment = models.CharField(max_length=150, blank=True, verbose_name = "Enviroment")
    test_temperature = models.CharField(max_length=150, blank=True, verbose_name = "Temp")
    test_time = models.CharField(max_length=25, blank=True, verbose_name = "Time")
    test_additive = models.CharField(max_length=150, blank=True, verbose_name = "Additive")
    subculture_type = models.CharField(max_length=25, blank=True, verbose_name = "Subculture/Seeding")
    incubation_time = models.CharField(max_length=25, blank=True, verbose_name = "Incubation Time")

    organism_id = models.ForeignKey(Organism, null=True, blank=True, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id")
    cell_id = models.ForeignKey(Cell, null=True, blank=True, verbose_name = "Cell ID", on_delete=models.DO_NOTHING,
        db_column="cell_id", related_name="%(class)s_cellid")
    
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    #------------------------------------------------
    class Meta:
        app_label = 'dscreen'
        db_table = 'assay'
        ordering=['assay_id']
        indexes = [
            models.Index(name="ass_aid_idx", fields=['assay_id']),
            models.Index(name="ass_aty_idx", fields=['assay_type']),
            models.Index(name="ass_sid_idx", fields=['sum_assay_id']),
            GinIndex(name="ass_pnl_idx", fields=['assay_panel']),
        ]

    # #------------------------------------------------
    # def __str__(self) -> str:
    #     return f"{self.run_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.run_id} [{self.run_type}]"

    #------------------------------------------------
    @classmethod
    def get(cls,AssayID,verbose=0):
        try:
            retInstance = cls.objects.get(assay_id=AssayID)
        except:
            if verbose:
                print(f"[Assay_ID Not Found] {AssayID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,AssayID,verbose=0):
        return cls.objects.filter(assay_id=AssayID).exists()
    
#-------------------------------------------------------------------------------------------------
class AssayData_MIC(CmpBatchList_Base):
    """
    List of MIC Values
    """
#-------------------------------------------------------------------------------------------------
    from dplate.models import TestPlate
#    from dorganism.models import Organism, Organism_Batch

    LIST_VIEW_FIELDS = {
        # "run_id":"Run ID",
        # "run_type":"Run Type",
        # "assay_note":"Assay",
        # "run_status":"Status",
        # "run_project":"Project",
        # "run_name":"Name",
        # "run_date":"Run Date",
        # "run_conditions":"Conditions",
        # "run_issues":"Issues",
    }

    DICTIONARY_FIELDS = {
        'pub_status':'Pub_Status',
        'data_quality':'Data_Quality',
    }
    
    # Primary Contraint
    testplate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "TestPlate ID", on_delete=models.DO_NOTHING,
        db_column="testplate_id", related_name="%(class)s_testplateid")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "TestWell ID")

    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")
    ora_assay_id = models.CharField(max_length=25, blank=True, verbose_name = "ora Assay ID")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 

    # Plate Information ----------------------------------------------------------------------------
    # test_date = models.DateField(null=True, blank=True, verbose_name = "Date")
    # plate_size = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Size", on_delete=models.DO_NOTHING,
    #     db_column="plate_size", related_name="%(class)s_platesize")
    # plate_material = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Material", on_delete=models.DO_NOTHING,
    #     db_column="plate_material", related_name="%(class)s_material")
    # orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
    #     db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    # test_media = models.CharField(max_length=40, blank=True, verbose_name = "Media")
    # test_dye = models.CharField(max_length=40, blank=True, verbose_name = "Dye")
    # test_additive = models.CharField(max_length=80, blank=True, verbose_name = "Additive")
    # readout_type = models.CharField(max_length=25, blank=True, verbose_name = "Readout Type")

    mic = models.CharField(max_length=50, verbose_name = "MIC")
    mic_unit = models.CharField(max_length=20, verbose_name = "Unit")
    mic_skips = models.SmallIntegerField(default=0, blank=True, verbose_name = "Skips")

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    analysis = models.CharField(max_length=15, verbose_name = "Analysis")

    inhibit_max = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMax")
    inhibit_min = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMin")
    conc_max = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMax")
    conc_min = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMin")
    n_conc = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Conc")

    data_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Quality", on_delete=models.DO_NOTHING,
        db_column="data_quality", related_name="%(class)s_dataquality")
    data_comment = models.CharField(max_length=50, blank=True, verbose_name = "Data Comment")
    valid = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Valid")

    ref_mic = models.CharField(max_length=150, blank=True, verbose_name = "Ref MIC")
    ref_mic_chk = models.SmallIntegerField(default=-1, blank=True, verbose_name = "d(Dilution)")
    ic50 = models.CharField(max_length=50, blank=True, verbose_name = "IC50")
    ic50_unit = models.CharField(max_length=20, blank=True, verbose_name = "IC50 Unit")
    ic50_pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "IC50 pScore")
    ic50_quality = models.CharField(max_length=20, blank=True, verbose_name = "IC50 Quality")
    ic50_r2 = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "IC50 r2")
    ic50_slope = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "IC50 Slope")

    pub_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pub Status", on_delete=models.DO_NOTHING,
        db_column="pub_status", related_name="%(class)s_pub_statust")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")
    chk_migration = models.SmallIntegerField(default=-1, blank=False, verbose_name = "Check for migration")

    class Meta:
        app_label = 'dscreen'
        db_table = 'assaydata_mic'
        ordering=['assay_id','testplate_id','testwell_id']
        constraints = [
            models.UniqueConstraint(name='assmic_loc_cst', fields=['testplate_id', 'testwell_id'], )
        ]        
        indexes = [
            GinIndex(name="assmic_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="assmic_rid_idx",fields=['run_id']),
            models.Index(name="assmic_sid_idx",fields=['assay_id']),
            models.Index(name="assmic_asc_idx",fields=['act_score']),
            models.Index(name="assmic_ana_idx",fields=['analysis']),
        #    models.Index(name="assmic_rot_idx",fields=['readout_type']),
            models.Index(name="assmic_act_idx",fields=['act_type']),
            models.Index(name="assmic_psc_idx",fields=['pscore']),
            models.Index(name="assmic_val_idx",fields=['valid']),
            models.Index(name="assmic_dqy_idx",fields=['data_quality']),
            models.Index(name="assmic_dcm_idx",fields=['data_comment']),
            models.Index(name="assmic_chkm_idx",fields=['chk_migration']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,PlateID,WellID,verbose=0):
        try:
            retInstance = cls.objects.get(testplate_id=PlateID, testwell_id=WellID)
        except:
            if verbose:
                logger.warning(f"[AssayData Not Found] {PlateID} {WellID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an AssayData  instance if found by name
    @classmethod
    def exists(cls,PlateID,WellID):
        return cls.objects.filter(testplate_id=PlateID, testwell_id=WellID).exists()

    #------------------------------------------------  
    def conv_list_to_string(self):
        super().conv_list_to_string()
        self.mic_lst        = COMPOUND_SEP.join([str(x) for x in self.mic if x > 0])
        self.mic_unit_lst   = COMPOUND_SEP.join([str(x) for x in self.mic_unit if x > 0])

    #------------------------------------------------  
    def conv_string_to_lst(self):
        super().conv_string_to_list()
        self.mic = strList_to_List(self.mic_lst,sep=COMPOUND_SEP,size=4,fill="")
        self.mic_unit = strList_to_List(self.mic_unit_lst,sep=COMPOUND_SEP,size=4,fill="")

   #------------------------------------------------
    def set_actscores(self,verbose=0):
        self.act_score = ActScore_DR(self.mic,self.mic_unit,DMax=self.inhibit_max)
        self.pscore = pScore(self.mic,self.mic_unit,self.inhibit_max,MW=self.cmpbatch_id.full_mw,gtShift=3,drMax2=40)
    
    #------------------------------------------------
    # 
    @classmethod    
    def calc_doseresponse(cls,):
        retInst = cls()
        
        return(retInst)    

#-------------------------------------------------------------------------------------------------
class AssayData_CC50(CmpBatchList_Base):
    """
    List of Cytotoxicty Values
    """
#-------------------------------------------------------------------------------------------------
    from dplate.models import TestPlate

    LIST_VIEW_FIELDS = {
        # "run_id":"Run ID",
    }

    DICTIONARY_FIELDS = {
        'pub_status':'Pub_Status',
        'data_quality':'Data_Quality',
    }
        
    # Primary Contraint
    testplate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "TestPlate ID", on_delete=models.DO_NOTHING,
        db_column="testplate_id", related_name="%(class)s_testplateid")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "TestWell ID")

    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")
    ora_assay_id = models.CharField(max_length=25, blank=True, verbose_name = "ora Assay ID")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 

    cc50 = models.CharField(max_length=50, blank=True, verbose_name = "CC50")
    cc50_unit = models.CharField(max_length=20, blank=True, verbose_name = "CC50 Unit")
    cc50_pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "CC50 pScore")
    cc50_quality = models.CharField(max_length=20, blank=True, verbose_name = "CC50 Quality")
    cc50_r2 = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "CC50 r2")
    cc50_slope = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CC50 Slope")

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    analysis = models.CharField(max_length=15, verbose_name = "Analysis")

    inhibit_max = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMax")
    inhibit_min = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMin")
    conc_max = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMax")
    conc_min = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMin")
    n_conc = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Conc")

    data_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Quality", on_delete=models.DO_NOTHING,
        db_column="data_quality", related_name="%(class)s_dataquality")
    data_comment = models.CharField(max_length=50, blank=True, verbose_name = "Data Comment")
    valid = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Valid")

    ref_cc50 = models.CharField(max_length=150, blank=True, verbose_name = "Ref MIC")
    ref_cc50_chk = models.SmallIntegerField(default=-1, blank=True, verbose_name = "d(Dilution)")

    pub_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pub Status", on_delete=models.DO_NOTHING,
        db_column="pub_status", related_name="%(class)s_pub_statust")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")
    chk_migration = models.SmallIntegerField(default=-1, blank=False, verbose_name = "Check for migration")

    class Meta:
        app_label = 'dscreen'
        db_table = 'assaydata_cc50'
        ordering=['assay_id','testplate_id','testwell_id']
        constraints = [
            models.UniqueConstraint(name='asscc50_loc_cst', fields=['testplate_id', 'testwell_id'], )
        ]        
        indexes = [
            GinIndex(name="asscc50_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="asscc50_rid_idx",fields=['run_id']),
            models.Index(name="asscc50_sid_idx",fields=['assay_id']),
            models.Index(name="asscc50_asc_idx",fields=['act_score']),
            models.Index(name="asscc50_ana_idx",fields=['analysis']),
        #    models.Index(name="asscc50_rot_idx",fields=['readout_type']),
            models.Index(name="asscc50_act_idx",fields=['act_type']),
            models.Index(name="asscc50_psc_idx",fields=['pscore']),
            models.Index(name="asscc50_val_idx",fields=['valid']),
            models.Index(name="asscc50_dqy_idx",fields=['data_quality']),
            models.Index(name="asscc50_dcm_idx",fields=['data_comment']),
            models.Index(name="asscc50_chkm_idx",fields=['chk_migration']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,PlateID,WellID,verbose=0):
        try:
            retInstance = cls.objects.get(testplate_id=PlateID, testwell_id=WellID)
        except:
            if verbose:
                logger.warning(f"[AssayData Not Found] {PlateID} {WellID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,PlateID,WellID):
        return cls.objects.filter(testplate_id=PlateID, testwell_id=WellID).exists()

    #------------------------------------------------  
    # def conv_list_to_string(self):
    #     super().conv_list_to_string()
    #     self.mic_lst        = COMPOUND_SEP.join([str(x) for x in self.mic if x > 0])
    #     self.mic_unit_lst   = COMPOUND_SEP.join([str(x) for x in self.mic_unit if x > 0])

    #------------------------------------------------  
    # def conv_string_to_lst(self):
    #     super().conv_string_to_list()
    #     self.mic = strList_to_List(self.mic_lst,sep=COMPOUND_SEP,size=4,fill="")
    #     self.mic_unit = strList_to_List(self.mic_unit_lst,sep=COMPOUND_SEP,size=4,fill="")

   #------------------------------------------------
    def set_actscores(self,verbose=0):
        self.act_score = ActScore_DR(self.cc50,self.cc50_unit,DMax=self.inhibit_max)
        self.pscore = pScore(self.cc50,self.cc50_unit,self.inhibit_max,MW=self.cmpbatch_id.full_mw,gtShift=3,drMax2=40)

#-------------------------------------------------------------------------------------------------
class AssayData_HC50(CmpBatchList_Base):
    """
    List of Haemolysis Values
    """
#-------------------------------------------------------------------------------------------------
    from dplate.models import TestPlate

    LIST_VIEW_FIELDS = {
        # "run_id":"Run ID",
    }

    DICTIONARY_FIELDS = {
        'pub_status':'Pub_Status',
        'data_quality':'Data_Quality',
    }
        
    # Primary Contraint
    testplate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "TestPlate ID", on_delete=models.DO_NOTHING,
        db_column="testplate_id", related_name="%(class)s_testplateid")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "TestWell ID")

    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")
    ora_assay_id = models.CharField(max_length=25, blank=True, verbose_name = "ora Assay ID")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 

    hc50 = models.CharField(max_length=50, blank=True, verbose_name = "HC50")
    hc50_unit = models.CharField(max_length=20, blank=True, verbose_name = "HC50 Unit")
    hc50_pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "HC50 pScore")
    hc50_quality = models.CharField(max_length=20, blank=True, verbose_name = "HC50 Quality")
    hc50_r2 = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "HC50 r2")
    hc50_slope = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "HC50 Slope")

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    hc10 = models.CharField(max_length=50, blank=True, verbose_name = "HC10")
    tox_type = models.CharField(max_length=5, blank=True, verbose_name = "Tox Type")
    tox_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Tox Score")

    analysis = models.CharField(max_length=15, verbose_name = "Analysis")

    inhibit_max = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMax")
    inhibit_min = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMin")
    conc_max = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMax")
    conc_min = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMin")
    n_conc = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Conc")

    data_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Quality", on_delete=models.DO_NOTHING,
        db_column="data_quality", related_name="%(class)s_dataquality")

#    data_quality = models.CharField(max_length=50, verbose_name = "Data Quality")
    data_comment = models.CharField(max_length=50, blank=True, verbose_name = "Data Comment")
    valid = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Valid")

    ref_hc50 = models.CharField(max_length=150, blank=True, verbose_name = "Ref MIC")
    ref_hc50_chk = models.SmallIntegerField(default=-1, blank=True, verbose_name = "d(Dilution)")

    pub_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pub Status", on_delete=models.DO_NOTHING,
        db_column="pub_status", related_name="%(class)s_pub_statust")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")
    chk_migration = models.SmallIntegerField(default=-1, blank=False, verbose_name = "Check for migration")

    class Meta:
        app_label = 'dscreen'
        db_table = 'assaydata_hc50'
        ordering=['assay_id','testplate_id','testwell_id']
        constraints = [
            models.UniqueConstraint(name='asshc50_loc_cst', fields=['testplate_id', 'testwell_id'], )
        ]        
        indexes = [
            GinIndex(name="asshc50_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="asshc50_rid_idx",fields=['run_id']),
            models.Index(name="asshc50_sid_idx",fields=['assay_id']),
            models.Index(name="asshc50_asc_idx",fields=['act_score']),
            models.Index(name="asshc50_tsc_idx",fields=['tox_score']),
            models.Index(name="asshc50_ana_idx",fields=['analysis']),
        #    models.Index(name="asshc50_rot_idx",fields=['readout_type']),
            models.Index(name="asshc50_act_idx",fields=['act_type']),
            models.Index(name="asshc50_tox_idx",fields=['tox_type']),
            models.Index(name="asshc50_psc_idx",fields=['pscore']),
            models.Index(name="asshc50_val_idx",fields=['valid']),
            models.Index(name="asshc50_dqy_idx",fields=['data_quality']),
            models.Index(name="asshc50_dcm_idx",fields=['data_comment']),
            models.Index(name="asshc50_chkm_idx",fields=['chk_migration']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,PlateID,WellID,verbose=0):
        try:
            retInstance = cls.objects.get(testplate_id=PlateID, testwell_id=WellID)
        except:
            if verbose:
                logger.warning(f"[AssayData Not Found] {PlateID} {WellID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,PlateID,WellID):
        return cls.objects.filter(testplate_id=PlateID, testwell_id=WellID).exists()

    #------------------------------------------------  
    # def conv_list_to_string(self):
    #     super().conv_list_to_string()
    #     self.mic_lst        = COMPOUND_SEP.join([str(x) for x in self.mic if x > 0])
    #     self.mic_unit_lst   = COMPOUND_SEP.join([str(x) for x in self.mic_unit if x > 0])

    #------------------------------------------------  
    # def conv_string_to_lst(self):
    #     super().conv_string_to_list()
    #     self.mic = strList_to_List(self.mic_lst,sep=COMPOUND_SEP,size=4,fill="")
    #     self.mic_unit = strList_to_List(self.mic_unit_lst,sep=COMPOUND_SEP,size=4,fill="")

   #------------------------------------------------
    def set_actscores(self,verbose=0):
        self.act_score = ActScore_DR(self.hc50,self.hc50_unit,DMax=self.inhibit_max,cutoff_inhib=10)
        self.pscore = pScore(self.hc50,self.hc50_unit,self.inhibit_max,MW=self.cmpbatch_id.full_mw,gtShift=3,drMax2=40)

#

class AssayData_SynergyMIC(CmpBatchList_Base):
    """
    List of Synergy 
    """
#-------------------------------------------------------------------------------------------------
    from dplate.models import TestPlate

    LIST_VIEW_FIELDS = {
        # "run_id":"Run ID",
    }

    DICTIONARY_FIELDS = {
        'pub_status':'Pub_Status',
        'data_quality':'Data_Quality',
    }

   # Primary Contraint
    testplate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "TestPlate ID", on_delete=models.DO_NOTHING,
        db_column="testplate_id", related_name="%(class)s_testplateid")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "TestWell ID")

    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 

    # Assay Data
    synmic = models.CharField(max_length=50, verbose_name = "MIC")
    synmic_unit = models.CharField(max_length=20, verbose_name = "Unit")
    synmic_skips = models.SmallIntegerField(default=0, blank=True, verbose_name = "Skips")

    fici = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "FICI")

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    analysis = models.CharField(max_length=15, verbose_name = "Analysis")

    inhibit_max = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMax")
    inhibit_min = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMin")
    # conc_max = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMax")
    # conc_min = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMin")
    n_conc = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Conc")

    # Data Quality
    data_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Quality", on_delete=models.DO_NOTHING,
        db_column="data_quality", related_name="%(class)s_dataquality")
    data_comment = models.CharField(max_length=50, blank=True, verbose_name = "Data Comment")
    valid = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Valid")

    # ref_mic = models.CharField(max_length=150, blank=True, verbose_name = "Ref MIC")
    # ref_mic_chk = models.SmallIntegerField(default=-1, blank=True, verbose_name = "d(Dilution)")
    # ic50 = models.CharField(max_length=50, blank=True, verbose_name = "IC50")
    # ic50_unit = models.CharField(max_length=20, blank=True, verbose_name = "IC50 Unit")
    # ic50_pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "IC50 pScore")
    # ic50_quality = models.CharField(max_length=20, blank=True, verbose_name = "IC50 Quality")
    # ic50_r2 = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "IC50 r2")
    # ic50_slope = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "IC50 Slope")

    pub_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pub Status", on_delete=models.DO_NOTHING,
        db_column="pub_status", related_name="%(class)s_pub_statust")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")


    class Meta:
        app_label = 'dscreen'
        db_table = 'assaydata_synmic'
        ordering=['assay_id','testplate_id','testwell_id']
        constraints = [
            models.UniqueConstraint(name='asssyn_loc_cst', fields=['testplate_id', 'testwell_id'], )
        ]        
        indexes = [
            GinIndex(name="asssyn_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="asssyn_rid_idx",fields=['run_id']),
            models.Index(name="asssyn_sid_idx",fields=['assay_id']),
            models.Index(name="asssyn_asc_idx",fields=['act_score']),
            models.Index(name="asssyn_ana_idx",fields=['analysis']),
        #    models.Index(name="asssyn_rot_idx",fields=['readout_type']),
            models.Index(name="asssyn_act_idx",fields=['act_type']),
            models.Index(name="asssyn_psc_idx",fields=['pscore']),
            models.Index(name="asssyn_val_idx",fields=['valid']),
            models.Index(name="asssyn_dqy_idx",fields=['data_quality']),
            models.Index(name="asssyn_dcm_idx",fields=['data_comment']),
        #    models.Index(name="asssyn_chkm_idx",fields=['chk_migration']),
        ]

#   &xID.                   Number,
#   Compounds               Varchar2(100),
##   TestPlate_ID            Varchar2(25),
##   TestWell_ID             Varchar2(3),
#   CompoundA1_ID           Varchar2(25),
#   CompoundA2_ID           Varchar2(25),
#   CompoundB1_ID           Varchar2(25),
#   CompoundB2_ID           Varchar2(25),
##   N_Compounds             Number(2,0),
##   AssayType_ID            Varchar2(25),
##   Test_Strain             Varchar2(25),
##   Test_Dye                Varchar2(25),
##   Test_Additive           Varchar2(25),
##   Test_Date               Date,
##   Run_ID                  Varchar2(25),
##   Analysis                Varchar2(10),
##   FICI_Value              Number(8,2),
#   FICI_Synergy            Number(3,0),
##   SYNMIC                  Varchar2(50),
##   SYNMIC_Unit             Varchar2(40),
##   SYNMIC_Value            Number,
##   SYNMIC_Prefix           Varchar2(2),
#   SYNMIC_CmpdA1           Number,
#   SYNMIC_CmpdA1_Unit      Varchar2(10),
#   SYNMIC_CmpdA2           Number,
#   SYNMIC_CmpdA2_Unit      Varchar2(10),
#   SYNMIC_CmpdB1           Number,
#   SYNMIC_CmpdB1_Unit      Varchar2(10),
#   SYNMIC_CmpdB2           Number,
#   SYNMIC_CmpdB2_Unit      Varchar2(10),
##   DMax                    Number(12,1),
##   DMin                    Number(12,1),
##   MIC_Skips               Number(3,0),
#   ConcA1_Min              Varchar2(20),
#   ConcA1_Max              Varchar2(20),
#   ConcB1_Min              Varchar2(20),
#   ConcB1_Max              Varchar2(20),
#   Checkerboard            Varchar2(10),
##   Data_Quality            Varchar2(20),
##   pScore                  Number(8,2),
#   MYSYC                   Number(3,0),
#   MYSYC_Synergy           Number(3,0),
#   MYSYC_Antagonism        Number(3,0),
#   MYSYC_beta              Number,
#   MYSYC_alpha12           Number,
#   MYSYC_alpha21           Number,
#   MYSYC_gamma12           Number,
#   MYSYC_gamma21           Number,
##   Hit                     Varchar2(4),
##   Active                  Varchar2(4),

# Assay (?)
# AssayData_SynMIC

class AssayData_SynergyIsobol(CmpBatchList_Base):
    """
    List of Synergy Isobolograms
    """
#-------------------------------------------------------------------------------------------------
    from dplate.models import TestPlate

    LIST_VIEW_FIELDS = {
        # "run_id":"Run ID",
    }

    DICTIONARY_FIELDS = {
        'pub_status':'Pub_Status',
        'data_quality':'Data_Quality',
    }

   # Primary Contraint
    testplate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "TestPlate ID", on_delete=models.DO_NOTHING,
        db_column="testplate_id", related_name="%(class)s_testplateid")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "TestWell ID")

    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 

    # Assay Data
    synmic = models.CharField(max_length=50, verbose_name = "MIC")
    synmic_unit = models.CharField(max_length=20, verbose_name = "Unit")
    synmic_skips = models.SmallIntegerField(default=0, blank=True, verbose_name = "Skips")

    fici = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "FICI")

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")
    pscore = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "pScore")

    analysis = models.CharField(max_length=15, verbose_name = "Analysis")

    inhibit_max = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMax")
    inhibit_min = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "DMin")
    # conc_max = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMax")
    # conc_min = models.DecimalField(max_digits=12, decimal_places=4, verbose_name = "CMin")
    n_conc = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Conc")

    # Data Quality
    data_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Quality", on_delete=models.DO_NOTHING,
        db_column="data_quality", related_name="%(class)s_dataquality")
    data_comment = models.CharField(max_length=50, blank=True, verbose_name = "Data Comment")
    valid = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Valid")

    class Meta:
        app_label = 'dscreen'
        db_table = 'assaydata_syniso'
        ordering=['assay_id','testplate_id','testwell_id']
        constraints = [
            models.UniqueConstraint(name='assiso_loc_cst', fields=['testplate_id', 'testwell_id'], )
        ]        
        indexes = [
            GinIndex(name="assiso_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="assiso_rid_idx",fields=['run_id']),
            models.Index(name="assiso_sid_idx",fields=['assay_id']),
            models.Index(name="assiso_asc_idx",fields=['act_score']),
        #    models.Index(name="assiso_ana_idx",fields=['analysis']),
        #    models.Index(name="assiso_rot_idx",fields=['readout_type']),
        #    models.Index(name="assiso_act_idx",fields=['act_type']),
        #    models.Index(name="assiso_psc_idx",fields=['pscore']),
        #    models.Index(name="assiso_val_idx",fields=['valid']),
        #    models.Index(name="assiso_dqy_idx",fields=['data_quality']),
        #    models.Index(name="assiso_dcm_idx",fields=['data_comment']),
        #    models.Index(name="assiso_chkm_idx",fields=['chk_migration']),
        ]

# AssayData_SynIsobol
#  &xID.                   Number,
#   Compounds               Varchar2(100),
#   CompoundA1_ID           Varchar2(25),
#   CompoundA2_ID           Varchar2(25),
#   CompoundB1_ID           Varchar2(25),
#   CompoundB2_ID           Varchar2(25),
#   N_Compounds             Number(2,0),
#   TestPlate_ID            Varchar2(25),
#   TestWell_ID             Varchar2(5),
#   AssayType_ID            Varchar2(25),
#   Test_Strain             Varchar2(25),
#   Test_Dye                Varchar2(25),
#   Test_Additive           Varchar2(25),
#   Test_Date               Date,
#   Run_ID                  Varchar2(25),
#   Analysis                Varchar2(10),
#   FICI_Value              Number(8,2),
#   FICI_Synergy            Number(3,0),
#   SYNMIC                  Varchar2(50),
#   SYNMIC_Unit             Varchar2(40),
#   SYNMIC_Value            Number,
#   SYNMIC_Prefix           Varchar2(2),
#   SYNMIC_CmpdA1           Number,
#   SYNMIC_CmpdA1_Unit      Varchar2(10),
#   SYNMIC_CmpdA2           Number,
#   SYNMIC_CmpdA2_Unit      Varchar2(10),
#   SYNMIC_CmpdB1           Number,
#   SYNMIC_CmpdB1_Unit      Varchar2(10),
#   SYNMIC_CmpdB2           Number,
#   SYNMIC_CmpdB2_Unit      Varchar2(10),
#   DMax                    Number(12,1),
#   DMin                    Number(12,1),
#   MIC_Skips               Number(3,0),
#   SC50                    Varchar2(20),
#   SC50_Prefix             Varchar2(2),
#   SC50_Value              Number,
#   SC50_Unit               Varchar2(10),
#   SC50_fSlope             Number,
#   SC50_fXC50              Number,
#   SC50_fR2                Number,
#   SC50_pScore             Number(8,2),
#   SC50_Quality            Varchar2(20),
#   ConcA1_Min              Varchar2(20),
#   ConcA1_Max              Varchar2(20),
#   ConcB1_Min              Varchar2(20),
#   ConcB1_Max              Varchar2(20),
#   Checkerboard            Varchar2(10),
#   Data_Quality            Varchar2(10),
#   Status                  Number(5),
#   aCreatedBy		      Varchar2(20),
#   aCreatedDate            Date,
#   aModifiedBy		      Varchar2(20),
#   aModifiedDate           Date,
#