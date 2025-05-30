#
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify
from django.contrib.auth.models import AbstractUser
from django.contrib.postgres.indexes import GinIndex

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from applib.data.str_lists import strList_to_List
from dcollab.models import Collab_Group, Collab_User
from dchem.models import Chem_Structure
from adjcoadd.constants import *
#

CMPBATCH_SOURCES = Choices( ('COADD','COADD CmpBatch'),
                          ('ABASE','ResearchGrp CmpBatch'),
                          ('LIBRARY','Library CmpBatch'),
                        )
import logging
logger = logging.getLogger(__name__)


#=================================================================================================
class Project(AuditModel):
    """
    List of Projects
    """
#=================================================================================================
    HEADER_FIELDS = {
        "project_id":{'Project ID': {'project_id':LinkList['project_id']}},
        "group_id.group_code":"Group",
        "group_id.country.name":"Country",
        "project_type":"Type",
        "project_status":"Status",
        "project_name":"Project Name",
        #"group_id":"Group",
        # "group_id.group_code":"Group",
    }

    DICTIONARY_FIELDS = {
        'project_type':'Project_Type',
        'project_status':'Project_Status',
        'provided_container':'Container_Type',
        'stock_conc_unit':'Unit_Concentration',
        'pub_status':'Pub_Status',
    }
    
    ID_SEQUENCE = 'Project'
    ID_PREFIX = 'P'
    ID_PAD = 5

    VIEW_GROUPS = [
        ['project_type','project_name','process_status','project_status','project_comment','received','completed',],
        ['provided_container','provided_comment','stock_status','stock_comment','stock_container','stock_conc','stock_conc_unit', ],
        ['compound_status','compound_comment','screen_status','screen_comment','data_status','data_comment',],
        ['report_status','report_comment','pub_status','pub_date','pub_name','source','source_code','reference']
    ]

    # Add Project Upload File Name

    project_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Project ID")
    project_name = models.CharField(max_length=150, blank=True, verbose_name = "Project Name")
    project_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Project Type", on_delete=models.DO_NOTHING,
        db_column="project_type", related_name="%(class)s_project_type")
    cpoz_id = models.CharField(max_length=50, blank=True, verbose_name = "CpOz ID")
    
    process_status = models.CharField(max_length=250, blank=True, verbose_name = "Process")
    project_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Project Status", on_delete=models.DO_NOTHING,
        db_column="project_status", related_name="%(class)s_project_status")
    project_comment = models.CharField(max_length=250, blank=True, verbose_name = "Project Comment")
    
    provided_container = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Provided Container", on_delete=models.DO_NOTHING,
        db_column="provided_container", related_name="%(class)s_provided_container")
    provided_comment = models.CharField(max_length=250, blank=True, verbose_name = "Provided Comment")
    
    received = models.DateField(null=True, blank=True, verbose_name = "Received")
    completed = models.DateField(null=True, blank=True, verbose_name = "Completed")
    
    stock_container = models.CharField(max_length=120, blank=True, verbose_name = "Stock Container")
    stock_conc = models.DecimalField(max_digits=9, decimal_places=2, default=0,verbose_name = "Stock Conc")
    stock_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Stock Conc Unit", on_delete=models.DO_NOTHING,
        db_column="stock_conc_unit", related_name="%(class)s_stock_conc_unit")
    stock_comment = models.CharField(max_length=150, blank=True, verbose_name = "Stock Comment")
    stock_status = ArrayField(models.CharField(max_length=50, null=True, blank=True), 
                                 size=20, verbose_name = "Stock Status", null=True, blank=True)
    
    compound_comment = models.CharField(max_length=150, blank=True, verbose_name = "Compound Comment")
    # Make Single Dictionary
    compound_status = ArrayField(models.CharField(max_length=50, null=True, blank=True), 
                                 size=20, verbose_name = "Compound Status", null=True, blank=True)
    
    screen_comment = models.CharField(max_length=150, blank=True, verbose_name = "Screen Comment")
    screen_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Screen Status", null=True, blank=True)

    # Make Single Dictionary
    data_comment = models.CharField(max_length=150, blank=True, verbose_name = "Data Comment")
    data_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Data Status", null=True, blank=True)

    report_comment = models.CharField(max_length=150, blank=True, verbose_name = "Report Comment")
    report_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Report Status", null=True, blank=True)

    group_id =  models.ForeignKey(Collab_Group, null=True, blank=True, verbose_name = "Project Owner", on_delete=models.DO_NOTHING,
        db_column="group_id", related_name="%(class)s_group_id")
    project_users =  ArrayField(models.CharField(max_length=25, null=True, blank=True), size=10, 
                             verbose_name = "Project Contacts", null=True, blank=True)
    #owner_users = models.ManyToManyField(Collab_User)
         
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    pub_name = models.CharField(max_length=150, blank=True, verbose_name = "Public Name")
    oldpub_status = models.CharField(max_length=200, null=True, blank=True, verbose_name = "Old Pub Status")
    # oldpub_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
    #                               size=20, verbose_name = "Public Status", null=True, blank=True)
    pub_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pub Status", on_delete=models.DO_NOTHING,
        db_column="pub_status", related_name="%(class)s_pub_statust")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")

    # -- Oracle-CastDB data ------------------------------------------------------------
    ORACLE_FIELDS = ['ora_project_id','ora_group_id','ora_contact_ids','ora_organisation','ora_psreport_date','ora_hcreport_date','ora_hvreport_date']

    ora_project_id = models.CharField(max_length=15, unique=True, verbose_name = "Old Project ID")
    ora_group_id = models.CharField(max_length=10, blank=True, verbose_name = "Old GroupID")
    ora_contact_ids = ArrayField(models.CharField(max_length=10, null=True, blank=True), size=2, 
                                 verbose_name = "Old ContactsUser", null=True, blank=True)
    ora_organisation = models.CharField(max_length=100, blank=True, verbose_name = "Old Organisation")
    ora_psreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="PS Report")
    ora_hcreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="HC Report")
    ora_hvreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="HV Report")


    # -- Calculated Fields - No View/Update ------------------------------------------
    CALCULATED_FIELDS = ['n_compounds','n_mcc_compounds','n_structure','n_barcode',
                         'n_motherplates','n_testplates','n_runids','n_assays',
                         'n_ps_compounds','n_dr_compounds','n_syn_compounds',
                         'n_sc_hits','n_mic_hits','n_tox_hits',
                         'screen_date']

    n_compounds = models.IntegerField(default=0, verbose_name = "#Cpmds")
    n_mcc_compounds = models.IntegerField(default=0, verbose_name = "#MCC")
    n_structure = models.IntegerField(default=0, verbose_name = "#Struc")
    n_barcode = models.IntegerField(default=0, verbose_name = "#BCode")
    n_motherplates = models.IntegerField(default=0, verbose_name = "#MP")
    n_testplates = models.IntegerField(default=0, verbose_name = "#TP")
    n_runids = models.IntegerField(default=0, verbose_name = "#Runs")
    n_assays = models.IntegerField(default=0, verbose_name = "#Assays")
    n_ps_compounds = models.IntegerField(default=0, verbose_name = "#PS")
    n_dr_compounds = models.IntegerField(default=0, verbose_name = "#DR")
    n_syn_compounds = models.IntegerField(default=0, verbose_name = "#SYN")
    n_sc_hits = models.IntegerField(default=0, verbose_name = "#Inhib Hits")
    n_mic_hits = models.IntegerField(default=0, verbose_name = "#MIC Hits")
    n_tox_hits = models.IntegerField(default=0, verbose_name = "#Tox Hits")
    screen_date = models.DateField(null=True, blank=True, verbose_name="Screen Date")


    class Meta:
        app_label = 'dsample'
        db_table = 'project'
        ordering=['project_id']
        indexes = [
            models.Index(name="prj_pname_idx", fields=['project_name']),
            models.Index(name="prj_opid_idx", fields=['ora_project_id']),
            models.Index(name="prj_ncmp_idx", fields=['n_compounds']),
            models.Index(name="prj_nstr_idx", fields=['n_structure']),
            models.Index(name="prj_nsh_idx", fields=['n_sc_hits']),
            models.Index(name="prj_nmh_idx", fields=['n_mic_hits']),
            models.Index(name="prj_nth_idx", fields=['n_tox_hits']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.project_id}  {self.source}"

    #------------------------------------------------
    @classmethod
    def get(cls,ProjectID,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            retInstance = cls.objects.get(project_id=ProjectID)
        except:
            retInstance = None
            if verbose:
                logger.warning(f"[Project Not Found] {ProjectID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,ProjectID,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        retValue = cls.objects.filter(project_id=ProjectID).exists()
        return(retValue)


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.project_id:
            self.project_id = self.next_id()
            if self.project_id: 
                super(Project, self).save(*args, **kwargs)
        else:
            super(Project, self).save(*args, **kwargs) 


#=================================================================================================
class Library(AuditModel):
    """
    List of Chem Library 
    """
#=================================================================================================
    DICTIONARY_FIELDS = {
        'project_type':'Project_Type',
        'project_status':'Project_Status',
        'provided_container':'Container_Type',
        'stock_conc_unit':'Unit_Concentration',
    }
    
    ID_SEQUENCE = 'Library'
    ID_PREFIX = 'L'
    ID_PAD = 5

    library_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Library ID")
    library_name = models.CharField(max_length=150, blank=True, verbose_name = "Name")
    library_version = models.CharField(max_length=15, blank=True, verbose_name = "Version")

    library_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="library_class", related_name="%(class)s_library_class")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    class Meta:
        app_label = 'dsample'
        db_table = 'library'
        ordering=['library_id']
        indexes = [
            models.Index(name="lib_lname_idx", fields=['library_name']),
            models.Index(name="lib_lclass_idx", fields=['library_class']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.library_id}  {self.source}"

    #------------------------------------------------
    @classmethod
    def get(cls,LibraryID=None,LibraryName=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if LibraryID:
                retInstance = cls.objects.get(library_id=LibraryID)
            elif LibraryName:
                retInstance = cls.objects.get(library_name=LibraryName)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if LibraryID:
                    logger.warning(f"[Library Not Found] {LibraryID} ")
                elif LibraryName:
                    logger.warning(f"[Library Not Found] {LibraryName} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,LibraryID=None,LibraryName=None,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        if LibraryID:
            retValue = cls.objects.filter(library_id=LibraryID).exists()
        elif LibraryName:
            retValue = cls.objects.filter(library_name=LibraryName).exists()
        else:
            retValue = False
        return(retValue)

#-------------------------------------------------------------------------------------------------
class Compound_Batch(AuditModel):
    """
    List of Compound Batches 
    """
#-------------------------------------------------------------------------------------------------
    DICTIONARY_FIELDS = {
        'batch_type':'CmpBatch_Type',
    }

    ID_SEQUENCE = 'CmpBatch'
    ID_PREFIX = 'CB'
    ID_PAD = 9
    
    cmpbatch_id = models.CharField(max_length=15, primary_key=True, verbose_name = "CmpBatch ID")

    batch_id  = models.CharField(default= '00',max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    batch_source = models.CharField(max_length=25, choices=CMPBATCH_SOURCES, blank=False, verbose_name = "Batch Source")
    batch_code = models.CharField(max_length=150, blank=True, verbose_name = "Batch Code")
#    batch_name = models.CharField(max_length=250, blank=True, verbose_name = "Batch Name")
#    batch_desc = models.CharField(max_length=512, blank=True, verbose_name = "Batch Description")
    
#    batch_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Batch Type", on_delete=models.DO_NOTHING,
#    db_column="batch_type", related_name="%(class)s_batch_type")
    
    previous_ids = models.CharField(max_length=100, blank=True, verbose_name = "Previous IDs")
    # parent_structure_ids = ArrayField(models.CharField(max_length=15, null=True, blank=True), size=4, verbose_name = "Panel", 
    #                                   null=True, blank=True)
    
    structure_type = models.CharField(max_length=400, blank=True, verbose_name = "Type")
    structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")

    # compound_id = models.ForeignKey(Compound, null=True, blank=True, verbose_name = "Compound ID", on_delete=models.DO_NOTHING,
    #     db_column="compound_id", related_name="%(class)s_compoundid")
    salt_code = models.CharField(max_length=120, blank=True, verbose_name = "Salts")
    smiles_extra = models.CharField(max_length=256, blank=True, verbose_name = "Smiles Salt")
    mw_extra = models.DecimalField(default=0, max_digits=12, decimal_places=3, verbose_name = "MW Salt")
    
    full_mw = models.FloatField(default=0, blank=True, verbose_name = "Full MW")
    full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")
    
    class Meta:
        app_label = 'dsample'
        db_table = 'cmpbatch'
        ordering=['cmpbatch_id']
        indexes = [
            models.Index(name="cmpbatch_src_idx", fields=['batch_source']),
            models.Index(name="cmpbatch_code_idx", fields=['batch_code']),
            models.Index(name="cmpbatch_stype_idx", fields=['structure_type']),
            models.Index(name="cmpbatch_fmw_idx", fields=['full_mw']),
            models.Index(name="cmpbatch_salt_idx", fields=['salt_code']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.cmpbatch_id}  {self.batch_code}"

    #------------------------------------------------
    @classmethod
    def get(cls,CmpBatchID,verbose=0):
    # Returns an instance by cmpbatch_id
        try:
            retInstance = cls.objects.get(cmpbatch_id=CmpBatchID)
        except:
            retInstance = None
            if verbose:
                logger.warning(f"[CmpBatch Not Found] {CmpBatchID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CmpBatchID,verbose=0):
    # Returns if an instance exists by cmpbatch_id
        retValue = cls.objects.filter(cmpbatch_id=CmpBatchID).exists()
        return(retValue)


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.cmpbatch_id:
            self.cmpbatch_id = self.next_id()
            if self.cmpbatch_id: 
                super(Compound_Batch, self).save(*args, **kwargs)
        else:
            super(Compound_Batch, self).save(*args, **kwargs) 


#-------------------------------------------------------------------------------------------------
class COADD_Compound(AuditModel):
    """
    List of CO-ADD Compounds as per Registration
    """
#-------------------------------------------------------------------------------------------------
    DICTIONARY_FIELDS = {
        'compound_type':'Compound_Type',
        'compound_source':'Compound_Source',
        'reg_amount_unit': 'Unit_Amount',
        'reg_volume_unit':'Unit_Volume',
        'reg_conc_unit':'Unit_Concentration',
    #    'stock_volume_unit':'Unit_Volume',
    }

    ID_SEQUENCE = 'COADD_Compound'
    ID_PREFIX = 'C'
    ID_PAD = 9
    
    compound_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Compound ID")
    compound_code = models.CharField(max_length=120, blank=True, verbose_name = "Code")
    compound_name = models.CharField(max_length=120, blank=True, verbose_name = "Name")
    compound_desc = models.CharField(max_length=150, blank=True, verbose_name = "Comment")

    project_id = models.ForeignKey(Project, null=True, blank=True, verbose_name = "Project ID", on_delete=models.DO_NOTHING,
        db_column="project_id", related_name="%(class)s_project_id")

    #cmpbatch_id = models.CharField(max_length=15, null=True, blank=True, verbose_name = "CmpBatch ID")
    cmpbatch_id = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")

    compound_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="compound_type", related_name="%(class)s_compound_type")

    compound_source = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="compound_source", related_name="%(class)s_compound_source")

    # compound_subtypes = ArrayField(models.CharField(max_length=50, null=True, blank=True), 
    #                                size=10, verbose_name = "Subtypes", null=True, blank=True)

    ora_compound_id = models.CharField(max_length=15, blank=True, verbose_name = "Old Compound ID")
    ora_project_id = models.CharField(max_length=15, blank=True, verbose_name = "Old Project ID")
    ora_compound_type = models.CharField(max_length=150, blank=True, verbose_name = "Old Compound Type")

    # CO-ADD - Registration ------
    reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")
    reg_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Reg MW")
    reg_mf = models.CharField(max_length=100, blank=True, verbose_name = "Reg MF")
    reg_structure = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Structure")
    reg_amount = models.DecimalField(max_digits=9, decimal_places=2, default=0, verbose_name = "Reg Amount")
    reg_amount_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Reg Amount Unit", on_delete=models.DO_NOTHING,
        db_column="reg_amount_unit", related_name="%(class)s_reg_amount_unit")
    reg_volume = models.DecimalField(max_digits=9, decimal_places=2, default=0, verbose_name = "Reg Volume")
    reg_volume_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Reg Volume Unit", on_delete=models.DO_NOTHING,
        db_column="reg_volume_unit", related_name="%(class)s_reg_volume_unit")
    reg_conc = models.DecimalField(max_digits=9, decimal_places=2, default=0,verbose_name = "Reg Conc")
    reg_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Reg Conc Unit", on_delete=models.DO_NOTHING,
        db_column="reg_conc_unit", related_name="%(class)s_reg_conc_unit")
    reg_solvent = models.CharField(max_length=100, blank=True, verbose_name = "Reg Solvent")
    
    # CO-ADD - Stock 
    prep_date = models.DateField(null=True, blank=True, verbose_name="Prepared")
    # stock_volume = models.DecimalField(max_digits=9, decimal_places=2, default=0, verbose_name = "Stock Volume")
    # stock_volume_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Stock Volume Unit", on_delete=models.DO_NOTHING,
    #     db_column="stock_amount_unit", related_name="%(class)s_stock_amount_unit")
    
    # CO-ADD - Strcuture Curation 
    std_status = models.CharField(max_length=10, blank=True, verbose_name = "Std Status")
    std_process = models.CharField(max_length=120, blank=True, verbose_name = "Std Process")
    std_issues = models.CharField(max_length=250, blank=True, verbose_name = "Std Issues")
    std_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Std Smiles")
    std_nfrag = models.SmallIntegerField(default=0, verbose_name = "Std nFrag")
    std_salt = models.CharField(max_length=100, blank=True, verbose_name = "Std Salt")
    std_ion = models.CharField(max_length=100, blank=True, verbose_name = "Std Ion")
    std_solvent = models.CharField(max_length=100, blank=True, verbose_name = "Std Solvent")
    std_metal = models.CharField(max_length=100, blank=True, verbose_name = "Std Metal")
    # std_structure_type = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
    #                              size=20, verbose_name = "Std Structure Type", null=True, blank=True)
    std_structure_type = models.CharField(max_length=400, blank=True, verbose_name = "Std Type")

    std_smiles_extra = models.CharField(max_length=256, blank=True, verbose_name = "Std Smiles Extra")
    std_mw = models.DecimalField(default=0, max_digits=12, decimal_places=3, verbose_name = "Std MW")
    std_mw_extra = models.DecimalField(default=0, max_digits=12, decimal_places=3, verbose_name = "Std MW Extra")
    std_mf = models.CharField(max_length=100, blank=True, verbose_name = "Std Total MF")

    # CO-ADD - Link to External ID's
    cpoz_sn = models.CharField(max_length=25, blank=True, verbose_name = "CpOz SN")
    cpoz_id = models.CharField(max_length=25, blank=True, verbose_name = "CpOz Lib ID")
    coadd_id = models.CharField(max_length=25, blank=True, verbose_name = "CO-ADD ID")
    chembl_id = models.CharField(max_length=25, blank=True, verbose_name = "ChEMBL ID")
    spark_id = models.CharField(max_length=25, blank=True, verbose_name = "SPARK ID")

    # CO-ADD - Publication
    pub_status = models.CharField(max_length=10, blank=True, verbose_name = "Pub Status")
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")

    class Meta:
        app_label = 'dsample'
        db_table = 'coadd_compound'
        ordering=['compound_id']
        indexes = [
            models.Index(name="coadd_name_idx", fields=['compound_name']),
            models.Index(name="coadd_code_idx", fields=['compound_code']),
            models.Index(name="coadd_type_idx", fields=['compound_type']),
            models.Index(name="coadd_pid_idx", fields=['project_id']),
            models.Index(name="coadd_cbid_idx", fields=['cmpbatch_id']),
            models.Index(name="coadd_ocid_idx", fields=['ora_compound_id']),
            models.Index(name="coadd_opid_idx", fields=['ora_project_id']),
            models.Index(name="coadd_sstat_idx", fields=['std_status']),
            models.Index(name="coadd_snfrag_idx", fields=['std_nfrag']),
            models.Index(name="coadd_sstyp_idx", fields=['std_structure_type']),
            models.Index(name="coadd_ssalt_idx", fields=['std_salt']),
            models.Index(name="coadd_smetal_idx", fields=['std_metal']),
            models.Index(name="coadd_pst_idx", fields=['pub_status']),   
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.compound_id}  {self.compound_code}"

    #------------------------------------------------
    @classmethod
    def get(cls,CompoundID,verbose=0):
    # Returns an instance by compound_id
        try:
            retInstance = cls.objects.get(compound_id=CompoundID)
        except:
            retInstance = None
            if verbose:
                logger.warning(f"[Compound Not Found] {CompoundID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CompoundID,verbose=0):
    # Returns if an instance exists by compound_id
        retValue = cls.objects.filter(compound_id=CompoundID).exists()
        return(retValue)

    #------------------------------------------------
    def save_batch(self, *args, **kwargs):
        if not self.cmpbatch_id:
            _batch = Compound_Batch.get(self.compound_id)
            if not _batch:
                _batch = Compound_Batch()
            _batch.cmpbatch_id= self.compound_id 
            _batch.batch_code = self.compound_code 
            #_batch.batch_type = self.compound_type
            
            # mw,mf, salt and structure_id
            _batch.save()
            self.cmpbatch_id = _batch
                

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.compound_id:
            self.compound_id = self.next_id()
            if self.compound_id:
                self.save_batch() 
                super(COADD_Compound, self).save(*args, **kwargs)
        else:
            self.save_batch()
            #logger.info(f" [Save COADD Compound] {self} {self.std_status}")
            super(COADD_Compound, self).save(*args, **kwargs) 


class ABase_Compound(AuditModel):
    """
    List of Abase Compounds as per Registration
    """
#-------------------------------------------------------------------------------------------------
    DICTIONARY_FIELDS = {
    }

    ID_SEQUENCE = 'ABase_Compound'
    ID_PREFIX = 'MCC'
    ID_PAD = 6

    compound_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Compound ID")
    compound_code = models.CharField(max_length=50, blank=True, verbose_name = "Code")
    compound_name = models.CharField(max_length=250, blank=True, verbose_name = "Name")
    compound_desc = models.CharField(max_length=250, blank=True, verbose_name = "Comment")

    reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")
    reg_molfile = models.TextField(max_length=15, blank=True, verbose_name = "Reg Molfile")
    reg_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Reg MW")
    reg_mf = models.CharField(max_length=100, blank=True, verbose_name = "Reg MF")
    
    structure_type = models.CharField(max_length=400, blank=True, verbose_name = "Type")
    structure_metal = models.CharField(max_length=100, blank=True, verbose_name = "Std Metal")
    structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")

    class Meta:
        app_label = 'dsample'
        db_table = 'abase_compound'
        ordering=['compound_id']
        indexes = [
            models.Index(name="abcmp_name_idx", fields=['compound_name']),
            models.Index(name="abcmp_code_idx", fields=['compound_code']),
            models.Index(name="abcmp_sid_idx", fields=['structure_id']),
        ]

    # @classmethod
    # def new_ABase_Compound_ID(cls,OldABaseID,verbose=0):
    #     return(OldABaseID.replace('MCC_','MCC'))


#-------------------------------------------------------------------------------------------------
    """
    List of Abase CmmpBatches as per Registration
    """
class ABase_Compound_Batch(AuditModel):

    DICTIONARY_FIELDS = {
        'init_amount_unit':'Unit_Amount',
    }

    # cmpbatch_id = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
    #     db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")
    cmpbatch_id = models.OneToOneField(Compound_Batch, primary_key=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
                                        db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")
    compound_id = models.ForeignKey(ABase_Compound, null=True, blank=True, verbose_name = "Compound ID", on_delete=models.DO_NOTHING,
        db_column="compound_id", related_name="%(class)s_compound_id")

    library_id = models.CharField(max_length=20, blank=True, verbose_name = "Library ID")
    project_id = models.ForeignKey(Project, null=True, blank=True, verbose_name = "Project ID", on_delete=models.DO_NOTHING,
            db_column="project_id", related_name="%(class)s_project_id")

    full_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Full MW")
    full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")
    salt_code = models.CharField(max_length=50, blank=True, verbose_name = "Salt Code")
    salt_equivalents = models.DecimalField(max_digits=7, decimal_places=2, default=0, verbose_name = "Salt Eq")
    solvate_code = models.CharField(max_length=50, blank=True, verbose_name = "Solvate Code")
    solvate_equivalents = models.DecimalField(max_digits=7, decimal_places=2, default=0, verbose_name = "Solvate Eq")
    conv_factor = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Conv Factor")
    supplier = models.CharField(max_length=50, blank=True, verbose_name = "Supplier")
    supplier_code = models.CharField(max_length=50, blank=True, verbose_name = "Supplier Code")
    supplier_batch = models.CharField(max_length=50, blank=True, verbose_name = "Supplier Batch")
    # supplier_po
    date_recieved = models.DateField(null=True, blank=True, verbose_name = "Received")
    init_amount = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Init Amount")
    init_amount_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Init Amount Unit", on_delete=models.DO_NOTHING,
        db_column="init_amount_unit", related_name="%(class)s_init_amount_unit")

    chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
        db_column="chemist", related_name="%(class)s_chemist")
    labbook_no = models.CharField(max_length=10, blank=True, verbose_name = "LabBook")
    labbook_page = models.CharField(max_length=10, blank=True, verbose_name = "LabBook Page")
    labbook_page_line = models.CharField(max_length=10, blank=True, verbose_name = "LabBook Page Line")

    class Meta:
        app_label = 'dsample'
        db_table = 'abase_cmpbatch'
        ordering=['compound_id']
        indexes = [
            models.Index(name="abcmpb_cmpb_idx", fields=['cmpbatch_id']),
            models.Index(name="abcmpb_salt_idx", fields=['salt_code']),
            models.Index(name="abcmpb_solv_idx", fields=['solvate_code']),
            models.Index(name="abcmpb_fmw_idx", fields=['full_mw']),
            models.Index(name="abcmpb_sup_idx", fields=['supplier']),
            models.Index(name="abcmpb_lid_idx", fields=['library_id']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,CompoundBatch, verbose=0):
    # Returns an instance by compound_id
        try:
            retInstance = cls.objects.get(cmpbatch_id=CompoundBatch)
        except:
            retInstance = None
            if verbose:
                logger.warning(f"[ABase CompoundBatch Not Found] {CompoundBatch} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CompoundBatch, verbose=0):
    # Returns if an instance exists by compound_id
        retValue = cls.objects.filter(cmpbatch_id=CompoundBatch).exists()
        return(retValue)


#-------------------------------------------------------------------------------------------------
class Library_Compound(AuditModel):
    """
    List of Library Compounds
    """
#-------------------------------------------------------------------------------------------------
    DICTIONARY_FIELDS = {
        'compound_type':'Compound_Type',
        # 'reg_amount_unit': 'Unit_Amount',
        # 'reg_volume_unit':'Unit_Volume',
        # 'reg_conc_unit':'Unit_Concentration',
    #    'stock_volume_unit':'Unit_Volume',
    }

    ID_SEQUENCE = 'Library_Compound'
    ID_PREFIX = 'LC'
    ID_PAD = 9
    
    compound_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Compound ID")
    compound_code = models.CharField(max_length=50, blank=True, verbose_name = "Code")
    
    compound_name = models.CharField(max_length=250, blank=True, verbose_name = "Name")
    compound_desc = models.CharField(max_length=250, blank=True, verbose_name = "Comment")

    library_id = models.ForeignKey(Library, null=True, blank=True, verbose_name = "Library ID", on_delete=models.DO_NOTHING,
        db_column="library_id", related_name="%(class)s_library_id")

    # cmpbatch_id = models.CharField(max_length=15, null=True, blank=True, verbose_name = "CmpBatch ID")
    cmpbatch_id = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")

    compound_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="compound_type", related_name="%(class)s_compound_type")
    
    reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")
    std_status = models.CharField(max_length=10, blank=True, verbose_name = "Std Status")
    std_process = models.CharField(max_length=120, blank=True, verbose_name = "Std Process")

    class Meta:
        app_label = 'dsample'
        db_table = 'library_compound'
        ordering=['compound_id']
        indexes = [
            models.Index(name="lcmp_code_idx", fields=['compound_code']),
            models.Index(name="lcmp_type_idx", fields=['compound_type']),
            models.Index(name="lcmp_lid_idx", fields=['library_id']),
            models.Index(name="lcmp_cbid_idx", fields=['cmpbatch_id']),
            models.Index(name="lcmp_sstat_idx", fields=['std_status']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,CompoundID, CompoundCode=None, LibraryID=None, verbose=0):
    # Returns an instance by compound_id
        try:
            if CompoundID:
                retInstance = cls.objects.get(compound_id=CompoundID)
            elif CompoundCode and LibraryID:
                retInstance = cls.objects.get(compound_code=CompoundCode,library_id=LibraryID)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if CompoundID:
                    logger.warning(f"[Compound Not Found] {CompoundID} ")
                elif CompoundCode and LibraryID:
                    logger.warning(f"[Compound Not Found] {CompoundCode} in {LibraryID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CompoundID, CompoundCode=None, LibraryID=None, verbose=0):
    # Returns if an instance exists by compound_id
        if CompoundID:
            retValue = cls.objects.filter(compound_id=CompoundID).exists()
        elif CompoundCode and LibraryID:
            retValue = cls.objects.filter(compound_code=CompoundCode,library_id=LibraryID).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.compound_id:
            self.compound_id = self.next_id()
            if self.compound_id:
                #self.save_batch() 
                super(Library_Compound, self).save(*args, **kwargs)
        else:
            #self.save_batch()
            super(Library_Compound, self).save(*args, **kwargs) 

#=================================================================================================
class Convert_ProjectID(AuditModel):
    """
    List of oraProjectID -> djProjectID
    """
#=================================================================================================

    ora_project_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Old Project ID")
    project_id = models.CharField(max_length=15,null=False, verbose_name = "Project ID")
    project_name = models.CharField(max_length=50, blank=True, verbose_name = "Project Name")

    class Meta:
        app_label = 'dsample'
        db_table = 'convert_projectid'
        ordering=['ora_project_id']
        indexes = [
            models.Index(name="wprj_pid_idx", fields=['project_id']),
            models.Index(name="wprj_pname_idx", fields=['project_name']),
#            models.Index(name="wprj_opid_idx", fields=['old_project_id']),
        ]


    @classmethod
    def new_COADD_Project_ID(cls,OldProjectID,verbose=0):

        if 'P' in OldProjectID:
            try:
                _cno = int(OldProjectID[1:])
            except:
                _cno = 0
                return(Project.str_id(_cno))
            if _cno > 0 :
                _newID = Project.str_id(_cno)
                newEntry = cls()
                newEntry.ora_project_id = OldProjectID
                newEntry.project_id = _newID
                newEntry.save()
            return(newEntry)


#=================================================================================================
class Convert_CompoundID(AuditModel):
    """
    List of oraCompoundID -> djCompoundID
    """
#=================================================================================================

    ora_compound_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Old Compound ID")
    compound_id = models.CharField(max_length=15,null=False, verbose_name = "Compound ID")
    compound_code = models.CharField(max_length=120, blank=True, verbose_name = "Code")
    compound_name = models.CharField(max_length=120, blank=True, verbose_name = "Name")
    project_id = models.CharField(max_length=15,null=False, verbose_name = "Project ID")
    sample_type = models.CharField(max_length=10,null=False, verbose_name = "Project ID")

    class Meta:
        app_label = 'dsample'
        db_table = 'convert_compoundid'
        ordering=['ora_compound_id']
        indexes = [
            models.Index(name="wcmpd_cid_idx", fields=['compound_id']),
            models.Index(name="wcmpd_ccode_idx", fields=['compound_code']),
            models.Index(name="wcmpd_stype_idx", fields=['sample_type']),
        ]

    @classmethod
    def new_COADD_Compound_ID(cls,OldCompoundID,verbose=0):

        if 'C0' in OldCompoundID:
            _cno = int(OldCompoundID[1:])
            _newID = COADD_Compound.str_id(_cno)

            newEntry = cls()
            newEntry.ora_compound_id = OldCompoundID
            newEntry.compound_id = _newID
            newEntry.save()
            return(newEntry)
        elif 'CM' in OldCompoundID:
            _cno = int(OldCompoundID[2:])+400000
            _newID = COADD_Compound.str_id(_cno)

            newEntry = cls()
            newEntry.ora_compound_id = OldCompoundID
            newEntry.compound_id = _newID
            newEntry.save()
            return(newEntry)
        
#-------------------------------------------------------------------------------------------------
class CmpBatchList_Base(AuditModel):    
#-------------------------------------------------------------------------------------------------

    STRING_FIELDS = ['cmpbatches']
    MAX_CMPBATCHES = 4

    cmpbatches = ""
    cmpbatch_lst = ArrayField(models.CharField(max_length=15, default=""), 
                                 size=MAX_CMPBATCHES, null=True, blank=True, db_index=True, verbose_name = "CmpBatch List")
    n_cmpbatches=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches")

    cmpbatch_id = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")

    class Meta:
        abstract = True
        ordering=['cmpbatch_lst']
        # To include in Child Models
        indexes = [
            GinIndex(name="cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="ncmpb_idx", fields=['n_cmpbatches']),
            models.Index(name="cmpbatch_idx", fields=['cmpbatch_lst']),
        ]

    #------------------------------------------------  
    def conv_list_to_string(self):
        if self.cmpbatch_lst:
            _CmpLst = [str(x) for x in self.cmpbatch_lst if x != ""]
            self.cmpbatches   = COMPOUND_SEP.join(_CmpLst)
            self.n_cmpbatches = len(_CmpLst)
        else:
            self.cmpbatches   = ''
            self.n_cmpbatches = 0

    #------------------------------------------------  
    def conv_string_to_list(self):
        self.cmpbatch_lst  = strList_to_List(self.cmpbatches,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def check_cmpbatch_id(self):
        _missing =[]
        for cmpbatch_id in [x for x in self.cmpbatch_lst if x != ""]:
            if not Compound_Batch.exists(cmpbatch_id):
                _missing.append(cmpbatch_id)
        if len(_missing) > 0:
            return({'Error': f"Compound_Batch not found {', '.join(_missing)}"})
        else:
            return(None)

    #------------------------------------------------  
    def set_cmpbatch_id(self,CmpBatchLst=None):
        if CmpBatchLst:
            _CmpLst = [str(x) for x in CmpBatchLst if x != ""]
            setattr(self,'cmpbatch_lst',_CmpLst)
        else:
           _CmpLst = getattr(self,'cmpbatch_lst')
        if _CmpLst: 
            setattr(self,'n_cmpbatches',len(_CmpLst))
            if self.n_cmpbatches == 1:
                self.cmpbatch_id = Compound_Batch.get(_CmpLst[0])
        else:
            setattr(self,'n_cmpbatches',0)
            self.cmpbatch_id = None 

    #------------------------------------------------  
    def __str__(self) -> str:
        return f"{self.cmpbatches}" 

    #------------------------------------------------  
    def __repr__(self) -> str:
        return f"{self.cmpbatches} " 

    #------------------------------------------------
    # Returns an User instance if found by name
    # @classmethod
    # def exists(cls,CmpBatchLst):
    #     if isinstance(CmpBatchLst,str):
    #         CmpBatchLst = [CmpBatchLst]
    #     return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst)).exists()

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def get(cls,CmpBatchLst,verbose=0):
        if isinstance(CmpBatchLst,str):
            CmpBatchLst = [CmpBatchLst]
        _qry = cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst))
        _cnt = _qry.count()
        if _cnt == 1:
            return(_qry[0])
        elif _cnt == 0:
            if verbose:
                logger.warning(f"[CmpBatchList Not Found] {CmpBatchLst} ")
            return(None)
        elif _cnt > 1:
            if verbose:
                logger.warning(f"[CmpBatchList Multiple Exist] {CmpBatchLst} : {_cnt}")
            return(None)
        return(None)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,CmpBatchLst):
        if isinstance(CmpBatchLst,str):
            CmpBatchLst = [CmpBatchLst]
        return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst)).exists()
    
    #------------------------------------------------
    # Clear/Reset cmpbatch data
    def clear_cmpbatch_data(self):
        self.cmpbatches = ""
        self.cmpbatch_lst = []
        self.n_cmpbatches=0
        self.cmpbatch_id = None
        
#-------------------------------------------------------------------------------------------------
class Sample_Base(CmpBatchList_Base):    
#-------------------------------------------------------------------------------------------------

    STRING_FIELDS = CmpBatchList_Base.STRING_FIELDS + ['concs','conc_units','conc_types'] 

    DICTIONARY_FIELDS = {
        'conc_unit_lst':'Unit_Concentration',
        'conc_type_lst':'Concentration_Type',
    }

    # MAX_CMPBATCHES = 4

    # cmpbatches = ""
    # cmpbatch_lst = ArrayField(models.CharField(max_length=15, default=""), 
    #                              size=MAX_CMPBATCHES, null=True, blank=True, db_index = True, verbose_name = "CmpBatch List")
    # n_cmpbatches=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches")

    concs = ""
    conc_lst = ArrayField(models.DecimalField(max_digits=9, decimal_places=4, default=0), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "Conc List", null=True, blank=True)
    conc_units = ""
    conc_unit_lst = ArrayField(models.CharField(max_length=10, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "ConcUnit List", null=True, blank=True)
    conc_types = ""
    conc_type_lst = ArrayField(models.CharField(max_length=5, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "ConcType List", null=True, blank=True)

    class Meta:
        abstract = True
        # ordering=['cmpbatch_lst']
        # # To include in Child Models
        # indexes = [
            # models.Index(name="conc_idx", fields=['conc_lst']),
            # models.Index(name="concuni_idx", fields=['conc_unit_lst']),
        # ]

    #------------------------------------------------  
    def conv_list_to_string(self):
        super().conv_list_to_string()
        self.concs        = ''
        self.conc_units   = ''
        self.conc_types   = ''
        
        if self.conc_lst:
            self.concs        = COMPOUND_SEP.join([str(x) for x in self.conc_lst if x > 0])
        if self.conc_unit_lst:
            self.conc_units   = COMPOUND_SEP.join([str(x) for x in self.conc_unit_lst if x != ""])
        if self.conc_type_lst:
            self.conc_types   = COMPOUND_SEP.join([str(x) for x in self.conc_type_lst if x != ""])

    #------------------------------------------------  
    def conv_string_to_list(self):
        super().conv_string_to_list()
        #self.cmpbatch_lst  = strList_to_List(self.cmpbatches,sep=COMPOUND_SEP,size=4,fill="")
        self.conc_lst      = strList_to_List(self.concs,sep=COMPOUND_SEP,size=4,fill=0)
        self.conc_unit_lst = strList_to_List(self.conc_units,sep=COMPOUND_SEP,size=4,fill="")
        self.conc_type_lst = strList_to_List(self.conc_types,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def check_conc_unit_dictionary(self):
        _missing = []
        for conc_unit in [x for x in self.conc_unit_lst if x != ""]:
            if not Dictionary.exists(self.DICTIONARY_FIELDS['conc_unit_lst'],conc_unit):
                _missing.append(conc_unit)
        if len(_missing) > 0:
            return({'Error': f"Conc_Unit not found {', '.join(_missing)}"})
        else:
            return(None)

    # #------------------------------------------------  
    # def check_cmpbatch_id(self):
    #     _missing =[]
    #     for cmpbatch_id in [x for x in self.cmpbatch_lst if x != ""]:
    #         if not Compound_Batch.exists(cmpbatch_id):
    #             _missing.append(cmpbatch_id)
    #     if len(_missing) > 0:
    #         return({'Error': f"Compound_Batch not found {', '.join(_missing)}"})
    #     else:
    #         return(None)

    # #------------------------------------------------  
    # def __str__(self) -> str:
    #     return f"{self.cmpbatches}" 

    #------------------------------------------------  
    def __repr__(self) -> str:
        return f"{self.cmpbatches} {self.concs} {self.conc_units}" 

    #------------------------------------------------
    # Clear/Reset cmpbatch data
    def clear_cmpbatch_data(self):
        super().clear_cmpbatch_data()
        self.concs = ""
        self.conc_lst = []
        self.conc_units = ""
        self.conc_unit_lst = []
        self.conc_types = ""
        self.conc_type_lst = []

    # #------------------------------------------------
    # # Returns an User instance if found by name
    # @classmethod
    # def get(cls,CmpBatchLst):
    #     pass

    # #------------------------------------------------
    # # Returns an User instance if found by name
    # @classmethod
    # def exists(cls,CmpBatchLst):
    #     return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst).exists()

