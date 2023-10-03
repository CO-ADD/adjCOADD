import re
from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from django.core.validators import RegexValidator


from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Organism Application Model
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Taxonomy(AuditModel):
    """
    Based on the NCBI Taxonomy at https://www.ncbi.nlm.nih.gov/taxonomy
    with available link https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=<Tax_ID>
    NCBI has an API but with Search/Fetch procedures https://www.ncbi.nlm.nih.gov/books/NBK25501/ 
    NCBI Taxonomy specific fields (could be defined as models.TextChoice classes):
        Tax_Rank        subspecies,species,species group,genus,order,class,family,phylum,no rank,varietas
        Division        Rodents, Bacteria, Mammals, Plants and Fungi, Primates
        Division_Code   ROD, BCT, MAM, PLN, PRI
    """
#=================================================================================================
    Choice_Dictionary = {
        'org_class':'Organism_Class',
        'division':'Organism_Division',
    }

    HEADER_FIELDS = {
        'organism_name':{'Organism Name': {'urlname': LinkList['taxonomny']}},  
        'tax_rank':'Rank',
        'org_class':'Class',
        'division':'Division', 
        'tax_id':{'Tax-ID': {'tax_id':LinkList['tax_id']}},
        'code':'Code', 
        'lineage':'Lineage', 
    }

    organism_name = models.CharField(primary_key=True, unique=True, max_length=100, verbose_name = "Name")
    urlname = models.SlugField(max_length=100, verbose_name = "URLName")
    other_names = models.CharField(max_length=100, blank=True, verbose_name = "Other Names")
    code = models.CharField(max_length=15, blank=True, verbose_name = "Code")
    org_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="org_class", related_name="%(class)s_class")
    tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")
    parent_tax_id = models.IntegerField(default=0, verbose_name = "NCBI Parent Tax ID") 
    tax_rank = models.CharField(max_length=50,  blank=True, verbose_name = "Taxonomy Rank")
    division = models.ForeignKey(Dictionary, null=True, verbose_name = "Division", on_delete=models.DO_NOTHING, 
        db_column="division", related_name='%(class)s_division')
    lineage = ArrayField(models.CharField(max_length=60, blank=True),size=30, null=True, verbose_name = "Lineage")
    
    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'taxonomy'
        ordering=['organism_name']
        indexes = [
            models.Index(name="tax_orgclass_idx", fields=['org_class']),
            models.Index(name="tax_taxid_idx", fields=['tax_id']),
            models.Index(name="tax_div_idx", fields=['division']),
            models.Index(name="tax_rnk_idx", fields=['tax_rank']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_name}"

    #------------------------------------------------
    @classmethod
    def get(cls,OrgName,verbose=0):
        try:
            retInstance = cls.objects.get(organism_name=OrgName.strip())
        except:
            if verbose:
                print(f"[Taxonomy Not Found] {OrgName} ")
            retInstance = None
        return(retInstance)
    #------------------------------------------------
    @classmethod
    def exists(cls,OrgName,verbose=0):
        return cls.objects.filter(organism_name=OrgName.strip()).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        self.urlname = slugify(self.organism_name,allow_unicode=False)
        super(Taxonomy, self).save()

        
#=================================================================================================
class Organism(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """
#=================================================================================================
    HEADER_FIELDS = {
#        'organism_name':{"VerboseName":'Organism Name','Updatable':False}
        'organism_id':{'Organism ID': {'organism_id':LinkList['organism_id']}}, 
        'organism_name':'Organism Name',
        'strain_ids':'Strain IDs',
        'sero_clone': 'Clone',
        'strain_type':'Strain Type',
        'strain_panel':'Panel',
        'strain_notes':'Notes',
        'res_property':'Phenotype',  
        'gen_property':'Genotype', 
        'strain_origin':'Origin',
        'source':"Source",
        'source_code':"Source Code",
        'reference': "Reference",
        'tax_id':{'Tax-ID': {'tax_id':LinkList['tax_id']}},
    }

    CARDS_FIELDS= {
        "source" : "Source",
        "strain_type": "Type",
        "res_property": "Phenotype",
        "gen_property": "Genotype",
    }

    FORM_GROUPS={
       'Group1': ["strain_ids", "sero_clone", "strain_code", "strain_type", "strain_panel", "strain_origin", "strain_notes"],
       'Group2': ["strain_identification", 'res_property','gen_property','oxygen_pref', 'source', 'source_code','reference','tax_id'],
       'Group3': ['mta_status','mta_notes','mta_document','risk_group','pathogen_group','lab_restriction','biologist'],
       'Group4': ['collect_date', 'collect_region', 'collect_country', 'collect_site', 'collect_specie', 'collect_tissue', 'patient_diagnosis', 'patient']
    }

    Choice_Dictionary = {
        'risk_group':'Risk_Group',
        'pathogen_group':'Pathogen_Group',
        'oxygen_pref':'Oxygen_Preference',
        'mta_status':'License_Status',
        'strain_type':'Strain_Type',
        'strain_panel':'Strain_Panel',
        'organism_class':'Organism_Class',
        'lab_restriction':'Lab_Restriction',
    }

    organism_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Organism ID") 
    organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
        db_column="organism_name", related_name="%(class)s_organism_name")
    strain_ids = models.CharField(max_length=200, blank=True, verbose_name = "Strain IDs") 
    strain_code= models.CharField(max_length=30, blank=True, verbose_name = "Strain Code")
    strain_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    strain_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    strain_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Strain Notes")
    res_property= models.CharField(max_length=1024, blank=True, verbose_name = "Phenotype")
    gen_property= models.CharField(max_length=1024, blank=True, verbose_name = "Genotype")
    sero_clone= models.CharField(max_length=126, blank=True, verbose_name = "MLST/Serotype")
    strain_identification = models.CharField(max_length=512, blank=True, verbose_name = "Strain Identification")
    strain_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Strain")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")
#    growth_preference = models.CharField(max_length=250, blank=True, verbose_name = "Growth/Screen Preference")
#    sequence_link = models.CharField(max_length=500, blank=True, verbose_name = "Sequence Link")

    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_mta")
    mta_notes = models.CharField(max_length=150, blank=True, verbose_name = "MTA Notes")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    mta_notes = models.CharField(max_length=512, blank=True, verbose_name = "MTA Notes")

    # received_date = models.DateField(null=True, blank=True, verbose_name = "Recieved")
    # received_as = models.CharField(max_length=120, blank=True, verbose_name = "Recieved as")
    # prep_notes= models.CharField(max_length=250, blank=True, verbose_name = "Preparation Notes")
    collect_date = models.DateField(null=True, blank=True, verbose_name = "Collection Date")
    # collect_notes = models.CharField(max_length=120, blank=True, verbose_name = "Notes")
    collect_region = models.CharField(max_length=25, blank=True, verbose_name = "Region/City")
    collect_country = models.CharField(max_length=25, blank=True, verbose_name = "Country")
    collect_site = models.CharField(max_length=50, blank=True, verbose_name = "Site/Org")

    collect_specie = models.CharField(max_length=20, blank=True, verbose_name = "From Specie/Location")
    collect_tissue = models.CharField(max_length=120, blank=True, verbose_name = "From Tissue/Organ")
    patient_diagnosis = models.CharField(max_length=120, blank=True, verbose_name = "Patient Diagnosis")
    patient = models.CharField(max_length=20, blank=True, verbose_name = "Patient Info")

    # collect_gender = models.CharField(max_length=5, blank=True, verbose_name = "Gender")
    # collect_age = models.CharField(max_length=10, blank=True, verbose_name = "Age")

    risk_group = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Risk Group", on_delete=models.DO_NOTHING,
        db_column="risk_group", related_name="%(class)s_risk")
    pathogen_group = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pathogen", on_delete=models.DO_NOTHING,
        db_column="pathogen_group", related_name="%(class)s_pathogen")
    lab_restriction = models.ForeignKey(Dictionary,null=True, blank=True, verbose_name = "Lab", on_delete=models.DO_NOTHING,
        db_column="lab_restriction", related_name="%(class)s_lab")
    oxygen_pref = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Oxygen", on_delete=models.DO_NOTHING,
        db_column="oxygen_pref", related_name="%(class)s_oxygen")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    # assoc_images = models.ManyToManyField(Image,verbose_name = "Images", blank=True,
    #     db_table = "org_img", related_name="%(class)s_image")
    assoc_documents = models.ManyToManyField(Document,verbose_name = "Douments", blank=True,
        db_table = "org_doc", related_name="%(class)s_document")

    #
    # List of Charfields with ';' list of URL like entries
    #  '<ExtSite>;<Text>;<Value1>;<Value2>' like ['nctc;NCTC 13368;13368','atcc;ATCC 90112;90112','chembl;; 
    # external_links = ArrayField(models.CharField(max_length=1200, blank=True),size=30, null=True, verbose_name = "Cross Reference") 
    #
    # or individual fields
    #
    # atcc	= models.URLField(blank=True, verbose_name = "atcc")	
    # nctc	= models.CharField(blank=True, verbose_name = "NCTC")	
    # cdc     = models.CharField(blank=True, verbose_name = "CDC ARBank")	


    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'organism'
        #ordering=['organism_id']
        indexes = [
            models.Index(name="org_stid_idx", fields=['strain_ids']),
            models.Index(name="org_stcode_idx", fields=['strain_code']),
            models.Index(name="org_strainid_idx", fields=['strain_type']),
            models.Index(name="org_stpanel_idx", fields=['strain_panel']),
            models.Index(name="org_source_idx", fields=['source']),
            models.Index(name="org_taxid_idx", fields=['tax_id']),
            models.Index(name="org_riskgrp_idx", fields=['risk_group']),
            models.Index(name="org_pathgrp_idx", fields=['pathogen_group']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_id} ({self.strain_code})"

    #------------------------------------------------
    @classmethod
    def str_OrganismID(cls,OrganimClass,OrganismNo) -> str:
    #
    # Input:    OrganismClass GN, GP,...
    #           OrganismNo 
    # Output:   Oragnism_ID as string like GN_0001 
    #
        return(f"{OrganimClass}{ORGANSIM_SEP}{OrganismNo:04d}")

    #------------------------------------------------
    @classmethod
    def exists(cls,OrgID,verbose=0):
    # Returns if an instance exists by organism_id
        return cls.objects.filter(organism_id=OrgID).exists()

    #------------------------------------------------
    @classmethod
    def get(cls,OrgID,verbose=0):
    # Returns an instance by organism_id
        try:
            retInstance = cls.objects.get(organism_id=OrgID)
        except:
            if verbose:
                print(f"[OrgansimID Not Found] {OrgID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def find_Next_OrganismID(cls,OrganismClass,OrganismClassTypes = ORGANISM_CLASSES) -> str:
        if OrganismClass in OrganismClassTypes:
            Organism_IDSq=Sequence(OrganismClass)
            Organism_nextID = next(Organism_IDSq)
            Organism_strID = cls.str_OrganismID(OrganismClass,Organism_nextID)
            while cls.exists(Organism_strID):
                Organism_nextID = next(Organism_IDSq)
                Organism_strID = cls.str_OrganismID(OrganismClass,Organism_nextID)
            return(Organism_strID)    
        else:
            return(None)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.organism_id: 
            self.organism_id = self.find_Next_OrganismID(str(self.organism_name.org_class.dict_value))
            if self.organism_id: 
                super(Organism, self).save(*args, **kwargs)
        else:
            super(Organism, self).save(*args, **kwargs) 

#=================================================================================================
class Organism_Batch(AuditModel):
    """
    Organism/Isolate Batch Collection
    """
#=================================================================================================
    HEADER_FIELDS = {
        "batch_id":"Batch ID",
        "batch_notes":"Batch Notes",
        "qc_status":"QC_Status",
        "qc_record": "QC Record",
        "stock_date":"Stock Date",
        "stock_level":"Stock Levels",
        "biologist":"Biologist"
    }

    Choice_Dictionary = {
        'qc_status':'QC_Status',
    }
    
    FORM_GROUPS = {
       'Group1': ["batch_id", "batch_notes", "qc_status", "qc_record", "stock_date", "stock_level", "biologist" ]
       }
    #SEP = '_'

    alphanumeric = RegexValidator(r'^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')

    orgbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "OrgBatch ID")
    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[alphanumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")
    qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC status", on_delete=models.DO_NOTHING,
        db_column="qc_status", related_name="%(class)s_qc")
    qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")
    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date") 
    stock_level = models.CharField(max_length=20, blank=True, verbose_name = "Stock Levels") 
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")
    
    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'orgbatch'
        ordering=['orgbatch_id']
        indexes = [
            models.Index(name="orgbatch_orgbatch_idx",fields=['organism_id','batch_id']),
            models.Index(name="orgbatch_qc_idx",fields=['qc_status']),
            # models.Index(name="orgbatch_supp_idx",fields=['supplier']),
            models.Index(name="orgbatch_sdate_idx",fields=['stock_date']),
            models.Index(name="orgbatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.orgbatch_id}"

    #------------------------------------------------
    @classmethod
    # Formats BatchNo:int -> BatchID:str 
    def str_BatchID(self,BatchNo:int) -> str:
        return(f"{BatchNo:02d}")
    #------------------------------------------------
    @classmethod
    # Formats OrganismID:str,BatchID:str -> OrgBatchID:str
    def str_OrgBatchID(self,OrganismID:str,BatchID:str) -> str:
        return(f"{OrganismID}{ORGBATCH_SEP}{BatchID}")

    #------------------------------------------------
    def find_Next_BatchID(self, OrganismID:str, BatchID:str=None) -> str:
        # Check for given BatchID    
        if BatchID:
            # Clean up BatchID - remove non alphanumeric character and make uppercase
            BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

            # Clean up BatchID - reformat numbers
            if BatchID.isnumeric():
                BatchID = self.str_BatchID(int(BatchID))

            next_OrgBatch = self.str_OrgBatchID(OrganismID,BatchID)
            if ~self.exists(next_OrgBatch):
                return(BatchID)

        # Find new BatchID    
        next_BatchNo = 1
        next_OrgBatch = self.str_OrgBatchID(OrganismID,self.str_BatchID(next_BatchNo))
        while self.exists(next_OrgBatch):
            next_BatchNo = next_BatchNo + 1
            next_OrgBatch = self.str_OrgBatchID(OrganismID,self.str_BatchID(next_BatchNo))
        return(self.str_BatchID(next_BatchNo))    

    #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,verbose=0):
    # Returns an instance if found by orgbatch_id
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID)
        except:
            if verbose:
                print(f"[OrgBatch Not Found] {OrgBatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(orgbatch_id=OrgBatchID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.orgbatch_id: 
            # creates new OrgBatchID
            OrgID = self.organism_id.organism_id
            BatchID = self.find_Next_BatchID(OrgID,self.batch_id)
            if BatchID:
                self.batch_id = BatchID
                self.orgbatch_id = self.str_OrgBatchID(OrgID,BatchID)
                super(Organism_Batch,self).save(*args, **kwargs)
        else:
            # confirms Batch_ID from OrgBatchID
            self.batch_id = str(self.orgbatch_id).replace(str(self.organism_id.organism_id),"").split(ORGBATCH_SEP)[1]
            super(Organism_Batch,self).save(*args, **kwargs)
        
# ================================================================================================

#-------------------------------------------------------------------------------------------------
class OrgBatch_Image(AuditModel):
#-------------------------------------------------------------------------------------------------
    HEADER_FIELDS = {
        'orgbatch_id': 'OrgBatch ID',
        'image_name':'Name', 
        'image_file':'Image',  
        'image_type':'Type',  
        'image_desc':'Description',
        'image_source':'Source',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    image_name =models.CharField(max_length=120, unique=True, verbose_name = "Name")
    image_file= models.ImageField(upload_to='images/orgbatch', verbose_name = "Image")
    image_type = models.CharField(max_length=25, verbose_name = "Image Type")
    image_desc = models.CharField(max_length=140, blank=True, verbose_name = "Description")
    image_source = models.CharField(max_length=50, blank=True, verbose_name = "Source")

    class Meta:
        app_label = 'dorganism'
        db_table = 'orgbatch_image'
        ordering=['orgbatch_id','image_name']
        indexes = [
            models.Index(name="obimg_name_idx",fields=['image_name']),
            models.Index(name="obimg_scr_idx",fields=['image_source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return str(self.pk)

    def __repr__(self) -> str:
        return f"[{self.image_name}] {str(self.orgbatch_id)} ({self.image_type})"

    #------------------------------------------------
    @classmethod
    def get(cls,ImgName,verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            retInstance = cls.objects.get(image_name=ImgName)
        except:
            if verbose:
                print(f"[OrgBatch Image Not Found] {ImgName} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,ImgName,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(image_name=ImgName).exists()

  

#=================================================================================================
class OrgBatch_Stock(AuditModel):
    """
    Stock of Organism/Isolate Batches
    
    """
#=================================================================================================
    HEADER_FIELDS={
        "orgbatch_id.orgbatch_id":{'OrgBatch ID': {'orgbatch_id.organism_id.organism_id':LinkList["organism_id"]}},
        "orgbatch_id.organism_id.organism_name":"Organism",
        #"orgbatch_id.organism_id.organism_name":{'Organism ID': {'orgbatch_id.organism_id.organism_id':LinkList['organism_id']}},
        "stock_type":"Stock Type",
        "n_created":"#Created",
        "n_left":"#Left",
        "location_freezer": "Freezer",
        "location_rack": "Rack",
        "location_column": "Column",
        "location_slot": "Slot",
        "stock_date":"Stock Date",
        "stock_note":"Stock Note",
        "biologist":"Biologist",
    }

    Choice_Dictionary = {
        'stock_type':'Stock_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    stock_type = models.ForeignKey(Dictionary, null=False, blank=False, verbose_name = "Stock Type", on_delete=models.DO_NOTHING,
        db_column="stock_type", related_name="%(class)s_stock")
    n_created = models.IntegerField(default=0, verbose_name = "#Created")
    n_left = models.IntegerField(default=0, verbose_name = "#Left")
    stock_date = models.DateField(null=True, blank = True, verbose_name = "Stock Date")
    stock_note = models.CharField(max_length=10, blank=True, verbose_name = "Stock Note")
    # passage_notes = models.CharField(max_length=30, blank=True, verbose_name = "Passage Notes")
    location_freezer = models.CharField(max_length=80, blank=True, verbose_name = "Freezer")
    location_rack = models.CharField(max_length=10, blank=True, verbose_name = "Rack")
    location_column = models.CharField(max_length=10, blank=True, verbose_name = "Column")
    location_slot = models.CharField(max_length=10, blank=True, verbose_name = "Slot")
    #stock_id = models.CharField(max_length=15, blank=True, verbose_name = "Stock ID")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'orgbatch_stock'
        ordering=['orgbatch_id','stock_type']
        indexes = [
            models.Index(name="orgbstock_stype_idx",fields=['stock_type']),
            models.Index(name="orgbstock_freezer_idx",fields=['location_freezer']),
            models.Index(name="orgbstock_stdate_idx",fields=['stock_date']),
            models.Index(name="orgbstock_nleft_idx",fields=['n_left']),

        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.pk} "
    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.orgbatch_id} {self.stock_type} {self.n_left}"

#    #------------------------------------------------
#     @classmethod
#     def get(cls,pkID,verbose=0):
#     # Returns an instance if found by orgbatch_id and stocktype
#         try:
#             retInstance = cls.objects.get(pk=pkID)
#         except:
#             if verbose:
#                 print(f"[OrgBatch Stock Not Found] {pkID}")
#             retInstance = None
#         return(retInstance)

#    #------------------------------------------------
#     @classmethod
#     def exists(cls,pkID,verbose=0):
#     # Returns if an instance exists by orgbatch_id and stocktype
#         return cls.objects.filter(pk=pkID).exists()


              
  
#=================================================================================================
class Organism_Culture(AuditModel):
    """
    Recommanded and optimised Growth/Culture conditions 
    
    """
#=================================================================================================
    HEADER_FIELDS = {
        # "organism_id":"Organism ID",
        "culture_type":"Type",
        "culture_source":"Source",
        "media":"Media",
        "addition":"Addition",
        "atmosphere":"Atmosphere",
        "temperature":"Temperature",
        "culture_notes":"Notes",
        "biologist":"Biologist"
    }

    Choice_Dictionary = {
        'culture_type':'Culture_Type',
        'culture_source':'Culture_Source',
    }
    
    FORM_GROUPS = {
        'Group1': ["culture_type", "culture_source", "media", "addition", "atmosphere", "temperature", "culture_notes", "biologist"]
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id")
    culture_type = models.ForeignKey(Dictionary, null=False, blank=False, verbose_name = "Culture Type", on_delete=models.DO_NOTHING,
        db_column="culture_type", related_name="%(class)s_culture_type")
    culture_source = models.ForeignKey(Dictionary, null=False, blank=False, verbose_name = "Source", on_delete=models.DO_NOTHING,
        db_column="culture_source", related_name="%(class)s_culture_source")
    media = models.CharField(max_length=120, blank=True, verbose_name = "Media") 
    addition = models.CharField(max_length=55, blank=True, verbose_name = "Addition") 
    atmosphere = models.CharField(max_length=120, blank=True, verbose_name = "Atmosphere") 
    temperature = models.CharField(max_length=25, blank=True, verbose_name = "Temperature") 
    # labware = models.CharField(max_length=120, blank=True, verbose_name = "Labware") 
    culture_notes = models.CharField(max_length=512,blank=True, verbose_name = "Notes") 
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'organism_culture'
        ordering=['organism_id','culture_type','media']
        indexes = [
            models.Index(name="orgcult_media_idx",fields=['media']),
            models.Index(name="orgcult_ctype_idx",fields=['culture_type']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.pk}"
    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.organism_id} {self.culture_type} {self.culture_source}"


    # # ------------------------------------------------


        
