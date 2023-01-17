from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from apputil.models import AuditModel, Dictionary, ApplicationUser
from psqlextra.indexes import UniqueIndex

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
#from django.utils.text import slugify
from apputil.utils import slugify

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Organism Application Model
#-------------------------------------------------------------------------------------------------

class Taxonomy(AuditModel):
    """
    Based on the NCBI Taxonomy at https://www.ncbi.nlm.nih.gov/taxonomy
    with available link https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=<Tax_ID>
    NCBI has an API but with Search/Fetch procedures https://www.ncbi.nlm.nih.gov/books/NBK25501/ 
    NCBI Taxonomy specific fields (could be defined as models.TextChoice classes):
        Tax_Rank        subspecies,species,species group,genus,order,class,family,phylum,no rank,varietas
        Division        Rodents, Bacteria, Mammals, Plants and Fungi, Primates
        Division_Code   ROD, BCT, MAM, PLN, PRI
    # """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'org_class':'Organism_Class',
        'division':'Organism_Division',
    }

    organism_name = models.CharField(primary_key=True, unique=True, max_length=100, verbose_name = "Specie")
    urlname = models.SlugField(max_length=100, verbose_name = "URLSpecie")
    other_names = models.CharField(max_length=100, blank=True, verbose_name = "Other Names")
    code = models.CharField(max_length=15, blank=True, verbose_name = "Code")
    org_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="org_class", related_name="Class+")
    tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")
    parent_tax_id = models.IntegerField(default=0, verbose_name = "NCBI Parent Tax ID") 
    tax_rank = models.CharField(max_length=50,  blank=True, verbose_name = "Taxonomy Rank")
    division = models.ForeignKey(Dictionary, null=True, verbose_name = "Division", on_delete=models.DO_NOTHING, 
        db_column="division", related_name='Division')
    lineage = ArrayField(models.CharField(max_length=60, blank=True),size=30, null=True)
    
    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'taxonomy'
        ordering=['organism_name']
        indexes = [
            models.Index(name="tax_orgclass_idx", fields=['org_class']),
            models.Index(name="tax_taxid_idx", fields=['tax_id']),
            models.Index(name="tax_div_idx", fields=['division']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_name}"

    #------------------------------------------------
    @classmethod
    def exists(self,OrgName,verbose=0):
        try:
            retInstance = self.objects.get(organism_name=OrgName.strip())
        except:
            if verbose:
                print(f"[Taxonomy Not Found] {OrgName} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        self.urlname = slugify(self.organism_name,lower=False,allow_unicode=False)
    #    self.urlname = slugify(self.organism_name,allow_unicode=False)
        super(Taxonomy, self).save()

    #------------------------------------------------
    def get_values(self, fields=TAXONOMY_FIELDs):
        value_list=super(Taxonomy, self).get_values(fields)
        return value_list
        
#-------------------------------------------------------------------------------------------------
class Organism(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """
#-------------------------------------------------------------------------------------------------
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

    #ORG_CLASSES = ['GN','GP','MB','FG','MA']
    #SEP = "_"

    organism_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Organism ID") 
    organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
        db_column="organism_name", related_name="organism_name+")
    strain_ids = models.CharField(max_length=200, blank=True, verbose_name = "Strain IDs")
    strain_code= models.CharField(max_length=30, blank=True, verbose_name = "Strain Code")
    strain_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    strain_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    res_property= models.CharField(max_length=350, blank=True, verbose_name = "Resistance Property")
    gen_property= models.CharField(max_length=350, blank=True, verbose_name = "Genetic Property")
    strain_origin = models.CharField(max_length=350, blank=True, verbose_name = "Origin of Strain")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")
    growth_preference = models.CharField(max_length=250, blank=True, verbose_name = "Growth/Screen Preference")
    strain_notes= models.CharField(max_length=250, blank=True, verbose_name = "Strain Notes")
    # tax_id contains a link, check django link field or class
    tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")
    sequence_link = models.CharField(max_length=500, blank=True, verbose_name = "Sequence Link")
    strain_identification = models.CharField(max_length=150, blank=True, verbose_name = "Strain Identification")
    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_MTA+")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    lab_restriction = models.ForeignKey(Dictionary,null=True, blank=True, verbose_name = "Lab", on_delete=models.DO_NOTHING,
        db_column="lab_restriction", related_name="%(class)s_MTA+")
    risk_group = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Risk Group", on_delete=models.DO_NOTHING,
        db_column="risk_group", related_name="%(class)s_Risk+")
    pathogen_group = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Pathogen", on_delete=models.DO_NOTHING,
        db_column="pathogen_group", related_name="%(class)s_Pathogen+")
    oxygen_pref = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Oxygen", on_delete=models.DO_NOTHING,
        db_column="oxygen_pref", related_name="%(class)s_Oxygen+")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_Biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'organism'
        ordering=['organism_id']
        indexes = [
            models.Index(name="org_stid_idx", fields=['strain_ids']),
            models.Index(name="org_stcode_idx", fields=['strain_code']),
            models.Index(name="org_strainid_idx", fields=['strain_type']),
            models.Index(name="org_stpanel_idx", fields=['strain_panel']),
            models.Index(name="org_taxid_idx", fields=['tax_id']),
            models.Index(name="org_riskgrp_idx", fields=['risk_group']),
            models.Index(name="org_pathgrp_idx", fields=['pathogen_group']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_id} ({self.strain_code})"

    #------------------------------------------------
    @classmethod
    def str_OrganismID(self,OrganimClass,OrganismNo) -> str:
    #
    # Input:    OrganismClass GN, GP,...
    #           OrganismNo 
    # Output:   Oragnism_ID as string like GN_0001 
    #
        return(f"{OrganimClass}{ORGANSIM_SEP}{OrganismNo:04d}")

    #------------------------------------------------
    @classmethod
    def exists(self,OrgID,verbose=0):
    #
    # Returns an instance by organism_id
    #
        try:
            retInstance = self.objects.get(organism_id=OrgID)
        except:
            if verbose:
                print(f"[OrgansimID Not Found] {OrgID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    def find_Next_OrganismID(self,OrganismClass,OrganismClassTypes = ORGANISM_CLASSES) -> str:
        if OrganismClass in OrganismClassTypes:
            Organism_IDSq=Sequence(OrganismClass)
            Organism_nextID = next(Organism_IDSq)
            Organism_strID = self.str_OrganismID(OrganismClass,Organism_nextID)
            while Organism.objects.filter(organism_id=Organism_strID).first():
                Organism_nextID = next(Organism_IDSq)
                Organism_strID = self.str_OrganismID(OrganismClass,Organism_nextID)
            return(Organism_strID)    
        else:
            return(None)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        #print(f"[save_Organism] ..")
        if not self.organism_id: #Object does not exists
            self.organism_id = self.find_Next_OrganismID(str(self.organism_name.org_class.dict_value))
            if self.organism_id: 
                super(Organism, self).save(*args, **kwargs)
        else:
            super(Organism, self).save(*args, **kwargs) 

    # ------------------------------------------------
    def get_values(self, fields=ORGANISM_FIELDs):
        value_list=super(Organism, self).get_values(fields)
        return value_list

#------------------------------------------------------------------------------------------------
class Organism_Batch(AuditModel):
    """
    Organism/Isolate Batch Collection
    """
#-------------------------------------------------------------------------------------------------

    Choice_Dictionary = {

        'qc_status':'QC_Status',
    }

    SEP = '_'

    orgbatch_id  = models.CharField(primary_key=True, max_length=10, verbose_name = "OrgBatch ID")
    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id+")
    batch_no  = models.IntegerField(default = -1, verbose_name = "Batch No")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")
    qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC Notes", on_delete=models.DO_NOTHING,
        db_column="qc_status", related_name="%(class)s_QC+")
    qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")
    supplier = models.CharField(max_length=250, blank=True, verbose_name = "Supplier")
    supplier_code = models.CharField(max_length=120, blank=True, verbose_name = "Supplier Code")
    supplier_po = models.CharField(max_length=120, blank=True, verbose_name = "Supplier PO")
    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date",editable=True) 
    stock_level = ArrayField(models.IntegerField(default=0), size=3, verbose_name = "Stock Levels", editable=False, default=list) 
    biologist = models.ForeignKey(ApplicationUser, null=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_Biologist")
    
    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'orgbatch'
        ordering=['orgbatch_id']
        indexes = [
            models.Index(name="orgbatch_orgbatch_idx",fields=['organism_id','batch_no']),
            models.Index(name="orgbatch_qc_idx",fields=['qc_status']),
            models.Index(name="orgbatch_supp_idx",fields=['supplier']),
            models.Index(name="orgbatch_sdate_idx",fields=['stock_date']),
            models.Index(name="orgbatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.orgbatch_id}"

    #------------------------------------------------
    def find_Next_BatchNo(self, OrganismID) -> int:
        next_BatchNo = 1
        while Organism_Batch.objects.filter(organism_id=OrganismID, batch_no=next_BatchNo).exists():
            next_BatchNo = next_BatchNo + 1
        return(next_BatchNo)    

    #------------------------------------------------
    @classmethod
    def str_OrgBatchID(self,OrganismID,BatchNo) -> str:
        return(f"{OrganismID}{ORGBATCH_SEP}{BatchNo:02d}")

    #------------------------------------------------
    @classmethod
    def exists(self,BatchID,verbose=0):
    #
    # Returns an instance if found by orgbatch_id
    #
        try:
            retInstance = self.objects.get(orgbatch_id=BatchID)
        except:
            if verbose:
                print(f"[OrgBatch Not Found] {BatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.orgbatch_id: #Object does not exists
            try:
                print(f'find orgbatchID with {self.organism_id.organism_id}')
            except Exception as err:
                print(err)
            Next_BatchNo = self.find_Next_BatchNo(self.organism_id.organism_id)
            if Next_BatchNo:
                self.batch_no = Next_BatchNo
                self.orgbatch_id = self.str_OrgBatchID(self.organism_id.organism_id,Next_BatchNo)
                super(Organism_Batch,self).save(*args, **kwargs)
        else:
            super(Organism_Batch,self).save(*args, **kwargs)
        
    # ------------------------------------------------
    def get_values(self, fields=ORGANISM_BATCH_FIELDs):
        value_list=super(Organism_Batch, self).get_values(fields)
        return value_list

#------------------------------------------------------------------------------------------------
class OrgBatch_Stock(AuditModel):
    """
    Stock of Organism/Isolate Batches
    
    """
#-------------------------------------------------------------------------------------------------

    Choice_Dictionary = {
        'stock_type':'Stock_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id+") 
    stock_note = models.CharField(max_length=10, blank=True, verbose_name = "Stock Note")
    passage_notes = models.CharField(max_length=30, blank=True, verbose_name = "Passage Notes")
    location_freezer = models.CharField(max_length=80, blank=True, verbose_name = "Freezer")
    location_rack = models.CharField(max_length=10, blank=True, verbose_name = "Rack")
    location_column = models.CharField(max_length=10, blank=True, verbose_name = "Column")
    location_slot = models.CharField(max_length=10, blank=True, verbose_name = "Slot")
    stock_type = models.ForeignKey(Dictionary, null=False, blank=False, verbose_name = "Stock Type", on_delete=models.DO_NOTHING,
        db_column="stock_type", related_name="%(class)s_stock_type+")
    stock_date = models.DateField(verbose_name = "Stock Date")
    stock_id = models.CharField(max_length=15, blank=True, verbose_name = "Stock ID")
    n_created = models.IntegerField(default=0, verbose_name = "#Vials created")
    n_left = models.IntegerField(default=0, verbose_name = "#Vials left")
    biologist = models.ForeignKey(ApplicationUser, null=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_Biologist")

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
            models.Index(name="orgbstock_stid_idx",fields=['stock_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.orgbatch_id} {self.stock_type} {self.n_left}"

   #------------------------------------------------
    @classmethod
    def exists(self,StockID,verbose=0):
    #
    # Returns an instance if found by orgbatch_id and stocktype
    #
        try:
            retInstance = self.objects.get(stock_id=StockID)
        except:
            if verbose:
                print(f"[OrgBatch Not Found] {StockID}")
            retInstance = None
        return(retInstance)

    # ------------------------------------------------
    def get_values(self, fields=ORGANISM_STOCK_FIELDs):
        value_list=super(OrgBatch_Stock, self).get_values(fields)
        return value_list

  
#-------------------------------------------------------------------------------------------------
class Organism_Culture(AuditModel):
    """
    Recommanded and optimised Growth/Culture conditions 
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'culture_type':'Culture_Type',
        'media_use':'Media_Use',
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id+")
    culture_type = models.CharField(max_length=120, blank=True, verbose_name = "Culture Type") 
    media_use = models.CharField(max_length=120, blank=True, verbose_name = "Media Use") 
    media = models.CharField(max_length=120, blank=True, verbose_name = "Media") 
    atmosphere = models.CharField(max_length=120, blank=True, verbose_name = "Atmosphere") 
    temperature = models.CharField(max_length=25, blank=True, verbose_name = "Temperature") 
    labware = models.CharField(max_length=120, blank=True, verbose_name = "Labware") 
    notes = models.CharField(max_length=512,blank=True, verbose_name = "Media") 
    biologist = models.ForeignKey(ApplicationUser, null=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_Biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dorganism'
        db_table = 'organism_culture'
        ordering=['organism_id','culture_type','media_use']
        indexes = [
            models.Index(name="orgcult_media_idx",fields=['media_use']),
            models.Index(name="orgcult_ctype_idx",fields=['culture_type']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_id} {self.media_use} {self.culture_type}"

    # ------------------------------------------------
    def get_values(self, fields=ORGANISM_CULTR_FIELDs):
        value_list=super(Organism_Culture, self).get_values(fields)
        return value_list
