from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from dcollab.models import Collab_Group, Collab_User
from dchem.models import Chem_Structure
from adjcoadd.constants import *



SAMPLE_SOURCES = Choices( ('COADD','COADD Sample'),
                          ('ABASE','ResearchGrp Sample'),
                          ('LIBRARY','Library Sample'),
                        )
#=================================================================================================
class Project(AuditModel):
    """
    List of Projects
    """
#=================================================================================================
    Choice_Dictionary = {
        'project_type':'Project_Type',
        'project_status':'Project_Status',
        'provided_container':'Container_Type',
        'stock_conc_unit':'Unit_Concentration',
    }
    
    ID_SEQUENCE = 'Project'
    ID_PREFIX = 'P'
    ID_PAD = 5

    project_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Project ID")
    project_name = models.CharField(max_length=150, blank=True, verbose_name = "Name")
    project_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="project_type", related_name="%(class)s_project_type")
    cpoz_id = models.CharField(max_length=50, blank=True, verbose_name = "CpOz ID")
    
    process_status = models.CharField(max_length=250, blank=True, verbose_name = "Process")
    project_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Project Status", on_delete=models.DO_NOTHING,
        db_column="project_status", related_name="%(class)s_project_status")
    project_comment = models.CharField(max_length=250, blank=True, verbose_name = "Comment")
    
    provided_container = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Provided Container", on_delete=models.DO_NOTHING,
        db_column="provided_container", related_name="%(class)s_provided_container")
    provided_comment = models.CharField(max_length=250, blank=True, verbose_name = "Comment")
    
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
    compound_status = ArrayField(models.CharField(max_length=50, null=True, blank=True), 
                                 size=20, verbose_name = "Compound Status", null=True, blank=True)
    
    screen_comment = models.CharField(max_length=150, blank=True, verbose_name = "Screen Comment")
    screen_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Screen Status", null=True, blank=True)

    data_comment = models.CharField(max_length=150, blank=True, verbose_name = "Data Comment")
    data_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Data Status", null=True, blank=True)

    report_comment = models.CharField(max_length=150, blank=True, verbose_name = "Report Comment")
    report_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Report Status", null=True, blank=True)

    owner_group =  models.ForeignKey(Collab_Group, null=True, blank=True, verbose_name = "Group", on_delete=models.DO_NOTHING,
        db_column="owner_group", related_name="%(class)s_owner_group")
    owner_user =  ArrayField(models.CharField(max_length=25, null=True, blank=True), size=10, 
                             verbose_name = "User", null=True, blank=True)
    #owner_users = models.ManyToManyField(Collab_User)

    ora_project_id = models.CharField(max_length=15, unique=True, verbose_name = "Old Project ID")
    ora_group_id = models.CharField(max_length=10, blank=True, verbose_name = "Old GroupID")
    ora_contact_ids = ArrayField(models.CharField(max_length=10, null=True, blank=True), size=2, 
                                 verbose_name = "Old ContactsUser", null=True, blank=True)
    ora_organisation = models.CharField(max_length=100, blank=True, verbose_name = "Old Organisation")
    ora_psreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="PS Report")
    ora_hcreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="HC Report")
    ora_hvreport_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="HV Report")
         
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    pub_name = models.CharField(max_length=150, blank=True, verbose_name = "Public Name")
    pub_status = ArrayField(models.CharField(max_length=20, null=True, blank=True), 
                                 size=20, verbose_name = "Public Status", null=True, blank=True)
    pub_date = models.DateField(null=True, blank=True,  editable=False, verbose_name="Published")

    class Meta:
        app_label = 'dsample'
        db_table = 'project'
        ordering=['project_id']
        indexes = [
            models.Index(name="prj_pname_idx", fields=['project_name']),
            models.Index(name="prj_opid_idx", fields=['ora_project_id']),
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
                print(f"[Project Not Found] {ProjectID} ")
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
    Choice_Dictionary = {
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
                    print(f"[Library Not Found] {LibraryID} ")
                elif LibraryName:
                    print(f"[Library Not Found] {LibraryName} ")
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
class Sample(AuditModel):
    """
    List of Samples and Compound Batches
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'sample_type':'Sample_Type',
    }

    ID_SEQUENCE = 'Sample'
    ID_PREFIX = 'S'
    ID_PAD = 9
    
    sample_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Sample ID")
    batch_id  = models.CharField(default= '00',max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    sample_source = models.CharField(max_length=25, choices=SAMPLE_SOURCES, blank=False, verbose_name = "Sample Source")
    sample_code = models.CharField(max_length=150, blank=True, verbose_name = "Sample Code")
#    sample_name = models.CharField(max_length=250, blank=True, verbose_name = "Sample Name")
#    sample_desc = models.CharField(max_length=512, blank=True, verbose_name = "Sample Description")
    
    sample_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Sample Type", on_delete=models.DO_NOTHING,
        db_column="sample_type", related_name="%(class)s_sample_type")
    
    previous_ids = models.CharField(max_length=100, blank=True, verbose_name = "Previous IDs")
    # parent_structure_ids = ArrayField(models.CharField(max_length=15, null=True, blank=True), size=4, verbose_name = "Panel", 
    #                                   null=True, blank=True)
    
    structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")

    # compound_id = models.ForeignKey(Compound, null=True, blank=True, verbose_name = "Compound ID", on_delete=models.DO_NOTHING,
    #     db_column="compound_id", related_name="%(class)s_compoundid")
    salt_code = models.CharField(max_length=120, blank=True, verbose_name = "Salts")
#     salt = models.CharField(max_length=500, blank=True, verbose_name = "Salt")
    full_mw = models.FloatField(default=0, blank=True, verbose_name = "Full MW")
    full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")
    
    class Meta:
        app_label = 'dsample'
        db_table = 'sample'
        ordering=['sample_id']
        indexes = [
            models.Index(name="sample_src_idx", fields=['sample_source']),
            models.Index(name="sample_code_idx", fields=['sample_code']),
            models.Index(name="sample_type_idx", fields=['sample_type']),
            models.Index(name="sample_fmw_idx", fields=['full_mw']),
            models.Index(name="sample_salt_idx", fields=['salt_code']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.sample_id}  {self.sample_code}"

    #------------------------------------------------
    @classmethod
    def get(cls,SampleID,verbose=0):
    # Returns an instance by sample_id
        try:
            retInstance = cls.objects.get(sample_id=SampleID)
        except:
            retInstance = None
            if verbose:
                print(f"[Sample Not Found] {SampleID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,SampleID,verbose=0):
    # Returns if an instance exists by sample_id
        retValue = cls.objects.filter(sample_id=SampleID).exists()
        return(retValue)


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.sample_id:
            self.sample_id = self.next_id()
            if self.sample_id: 
                super(Sample, self).save(*args, **kwargs)
        else:
            super(Sample, self).save(*args, **kwargs) 





#-------------------------------------------------------------------------------------------------
class COADD_Compound(AuditModel):
    """
    List of CO-ADD Samples as per Registration
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'compound_type':'Sample_Type',
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

    sample_id = models.ForeignKey(Sample, null=True, blank=True, verbose_name = "Sample ID", on_delete=models.DO_NOTHING,
        db_column="sample_id", related_name="%(class)s_sample_id")

    compound_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="compound_type", related_name="%(class)s_compound_type")

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
    std_action = models.CharField(max_length=120, blank=True, verbose_name = "Std Action")
    std_process = models.CharField(max_length=120, blank=True, verbose_name = "Std Process")
    std_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Std Smiles")
    std_salt = models.CharField(max_length=2048, blank=True, verbose_name = "Std Salt")
    std_ion = models.CharField(max_length=2048, blank=True, verbose_name = "Std Ion")
    std_solvent = models.CharField(max_length=2048, blank=True, verbose_name = "Std Solvent")

    std_mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Std MW")
    std_mf = models.CharField(max_length=100, blank=True, verbose_name = "Std MF")

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
            models.Index(name="coadd_sid_idx", fields=['sample_id']),
            models.Index(name="coadd_ocid_idx", fields=['ora_compound_id']),
            models.Index(name="coadd_opid_idx", fields=['ora_project_id']),
            models.Index(name="coadd_sst_idx", fields=['std_status']),
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
                print(f"[Sample Not Found] {CompoundID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CompoundID,verbose=0):
    # Returns if an instance exists by compound_id
        retValue = cls.objects.filter(compound_id=CompoundID).exists()
        return(retValue)

    #------------------------------------------------
    def save_sample(self, *args, **kwargs):
        if not self.sample_id:
            _sample = Sample.get(self.compound_id)
            if not _sample:
                _sample = Sample()
            _sample.sample_id = self.compound_id 
            _sample.sample_code = self.compound_code 
            _sample.sample_type = self.compound_type
            
            # mw,mf, salt and structure_id
            _sample.save()
            self.sample_id = _sample
                

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.compound_id:
            self.compound_id = self.next_id()
            if self.compound_id:
                self.save_sample() 
                super(COADD_Compound, self).save(*args, **kwargs)
        else:
            self.save_sample()
            super(COADD_Compound, self).save(*args, **kwargs) 



#-------------------------------------------------------------------------------------------------
class Library_Compound(AuditModel):
    """
    List of CO-ADD Samples as per Registration
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
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

    sample_id = models.ForeignKey(Sample, null=True, blank=True, verbose_name = "Sample ID", on_delete=models.DO_NOTHING,
        db_column="sample_id", related_name="%(class)s_sample_id")

    compound_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="compound_type", related_name="%(class)s_compound_type")
    
    reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")

    class Meta:
        app_label = 'dsample'
        db_table = 'library_compound'
        ordering=['compound_id']
        indexes = [
            models.Index(name="lcmp_code_idx", fields=['compound_code']),
            models.Index(name="lcmp_type_idx", fields=['compound_type']),
            models.Index(name="lcmp_lid_idx", fields=['library_id']),
            models.Index(name="lcmp_sid_idx", fields=['sample_id']),
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
                    print(f"[Compound Not Found] {CompoundID} ")
                elif CompoundCode and LibraryID:
                    print(f"[Compound Not Found] {CompoundCode} in {LibraryID} ")
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
                #self.save_sample() 
                super(Library_Compound, self).save(*args, **kwargs)
        else:
            #self.save_sample()
            super(Library_Compound, self).save(*args, **kwargs) 



# #=================================================================================================
# class Sample_Batch(AuditModel):
#     """
#     List of Sample Batches 
#     """
# #=================================================================================================
#     ID_SEQUENCE = 'Sample'
#     ID_PREFIX = 'SB'
#     ID_PAD = 9

#     samplebatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "SampleBatch ID")
#     sample_id = models.ForeignKey(Sample, null=True, blank=True, verbose_name = "Sample ID", on_delete=models.DO_NOTHING,
#         db_column="sample_id", related_name="%(class)s_sample_id")
#     previous_batch_id= models.CharField(max_length=20, blank=True, verbose_name = "Previous SampleBatch ID")
#     batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
#     batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    # structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
    #     db_column="structure_id", related_name="%(class)s_structureid")

#     salt = models.CharField(max_length=500, blank=True, verbose_name = "Salt")
#     mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
#     mw = models.FloatField(default=0, blank=True, verbose_name ="MW")

#     quality_source = models.CharField(max_length=150, blank=True, verbose_name = "QC Source")
#     qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC status", on_delete=models.DO_NOTHING,
#         db_column="qc_status", related_name="%(class)s_qc")
#     qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")
#     stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date") 
#     stock_level = models.CharField(max_length=20, blank=True, verbose_name = "Stock Levels") 
#     chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
#         db_column="chemist", related_name="%(class)s_chemist")

#     #------------------------------------------------
#     class Meta:
#         app_label = 'dsample'
#         db_table = 'samplebatch'
#         ordering=['samplebatch_id']
#         indexes = [
#             models.Index(name="smpbatch_samplebatch_idx",fields=['sample_id','batch_id']),
#             models.Index(name="smpbatch_qc_idx",fields=['qc_status']),
#             models.Index(name="smpbatch_sdate_idx",fields=['stock_date']),
#             models.Index(name="smpbatch_slevel_idx",fields=['stock_level']),
#         ]

#     #------------------------------------------------
#     def __str__(self) -> str:
#         return f"{self.samplebatch_id}"

#     #------------------------------------------------
#     @classmethod
#     # Formats BatchNo:int -> BatchID:str 
#     def str_BatchID(self,BatchNo:int) -> str:
#         return(f"{BatchNo:02d}")
#     #------------------------------------------------
#     @classmethod
#     def str_SampleBatchID(self,SampleID:str,BatchID:str) -> str:
#         return(f"{SampleID}{SAMPLEBATCH_SEP}{BatchID}")

#     #------------------------------------------------
#     def find_Next_BatchID(self, SampleID:str, BatchID:str=None) -> str:
#         # Check for given BatchID    
#         if BatchID:
#             # Clean up BatchID - remove non alphanumeric character and make uppercase
#             BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

#             # Clean up BatchID - reformat numbers
#             if BatchID.isnumeric():
#                 BatchID = self.str_BatchID(int(BatchID))

#             next_SampleBatch = self.str_SampleBatchID(SampleID,BatchID)
#             if not self.exists(next_SampleBatch):
#                 return(BatchID)

#         # Find new BatchID    
#         next_BatchNo = 1
#         next_SampleBatch = self.str_SampleBatchID(SampleID,self.str_BatchID(next_BatchNo))
#         while self.exists(next_SampleBatch):
#             next_BatchNo = next_BatchNo + 1
#             next_SampleBatch = self.str_SampleBatchID(SampleID,self.str_BatchID(next_BatchNo))
#         return(self.str_BatchID(next_BatchNo))    

#     #------------------------------------------------
#     @classmethod
#     def get(cls,SampleBatchID,verbose=0):
#     # Returns an instance if found by samplebatch_id
#         try:
#             retInstance = cls.objects.get(samplebatch_id=SampleBatchID)
#         except:
#             if verbose:
#                 print(f"[SampleBatch Not Found] {SampleBatchID} ")
#             retInstance = None
#         return(retInstance)

#     #------------------------------------------------
#     @classmethod
#     def exists(cls,SampleBatchID,verbose=0):
#     # Returns if instance exists
#         return cls.objects.filter(samplebatch_id=SampleBatchID).exists()

#     #------------------------------------------------
#     def save(self, *args, **kwargs):
#         if not self.samplebatch_id: 
#             # creates new SampleBatchID
#             SampleID = self.sample_id.sample_id
#             BatchID = self.find_Next_BatchID(SampleID,self.batch_id)
#             if BatchID:
#                 self.batch_id = BatchID
#                 self.samplebatch_id = self.str_SampleBatchID(SampleID,BatchID)
#                 super(Sample_Batch,self).save(*args, **kwargs)
#         else:
#             # confirms Batch_ID from SampleBatchID
#             self.batch_id = str(self.samplebatch_id).replace(str(self.sample_id.sample_id),"").split(SAMPLEBATCH_SEP)[1]
#             super(Sample_Batch,self).save(*args, **kwargs)
#             #print(f"[SampleBatch.save]: {self.samplebatch_id}")
    
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
