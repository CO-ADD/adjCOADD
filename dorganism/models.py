from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from apputil.models import AuditModel, Dictionary, ApplicationUser
from psqlextra.indexes import UniqueIndex

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError

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

    organism_name = models.CharField(primary_key=True, unique=True, max_length=100, verbose_name = "Specie")
    other_names = models.CharField(blank=True, max_length=100, verbose_name = "Other Names")
    code = models.CharField(blank=True, max_length=15, verbose_name = "Code")
    org_class = models.ForeignKey(Dictionary, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    tax_id = models.IntegerField(verbose_name = "NCBI Tax ID")
    parent_tax_id = models.IntegerField(verbose_name = "NCBI Parent Tax ID") #empty from no.207668 
    tax_rank = models.CharField(blank=True, max_length=50, verbose_name = "Taxonomy Rank")
    division = models.ForeignKey(Dictionary, blank=True, null=True, verbose_name = "Division", related_name='requests_Div', on_delete=models.DO_NOTHING)
    lineage = ArrayField(models.CharField(max_length=25, null=True, blank=True),size = 25)
    
    #------------------------------------------------
    class Meta:
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


#-------------------------------------------------------------------------------------------------
class Organism(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'risk_group':'Risk_Group',
        'pathogen_group':'Pathogen_Group',
        'bio_approval':'Biol_Approval',
        'oxygen_pref':'Oxygen_Preference',
        'mta_status':'License_Status',
        'strain_type':'Strain_Type',
        'strain_panel':'Strain_Panel',
        'organism_class':'Organism_Class',
    }


    organism_id = models.CharField(primary_key=True, unique=True, blank=True, max_length=15, verbose_name = "Organism ID") #be blank for automatic generate a new one?
    organism_name= models.ForeignKey(Taxonomy, verbose_name = "Organism Name", on_delete=models.DO_NOTHING) #models do nothing?
    #Organism_desc= models.CharField(blank=True, max_length=150, verbose_name = "Organism Description")
    strain_id = models.CharField(blank=False, null=False, max_length=50, db_index=True, verbose_name = "Strain ID")
    strain_otherids = models.CharField(blank=True, null=True, max_length=150, verbose_name = "Strain OtherID")
    strain_code= models.CharField(blank=True, null=True, max_length=15, verbose_name = "Strain Code")
    #strain_Desc= models.CharField(blank=True, max_length=150, verbose_name = "Strain Description")
    strain_notes= models.CharField(blank=True, null=True, max_length=250, verbose_name = "Strain Notes")
    strain_origin = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Origin")
    strain_property = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Property")
    strain_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    strain_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    strain_tissue= models.CharField(blank=True, null=True, max_length=250, verbose_name = "Strain Tissue")
    biologist = models.ForeignKey(ApplicationUser, verbose_name = "Biologist", on_delete=models.DO_NOTHING) #models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")
    sequence = models.CharField(blank=True, null=True, max_length=100, verbose_name = "Sequence")
    sequence_link = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Sequence Link")
    tax_id = models.IntegerField(null=True, verbose_name = "NCBI Tax ID", default=0)
    risk_group = models.CharField(blank=True, null=True, max_length=50,verbose_name = "Risk Group")
    pathogen_group = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Pathogen",)
    import_permit = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Import Permit")
    bio_approval = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Biological Approval") 
    special_precaution = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Special Precaution")
    lab_restriction = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Lab")
    mta_status = models.CharField(blank=True, null=True, max_length=150,verbose_name = "MTA Status") 
    mta_document = models.CharField(blank=True, null=True, max_length=150, verbose_name = "MTA Document")
    oxygen_pref = models.CharField(blank=True, null=True, max_length=250, verbose_name = "Oxygen")
    atmosphere_pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Atmosphere")
    nutrient_pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Nutirent")
    biofilm_pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Biofilm")

          #------------------------------------------------
    class Meta:
        db_table = 'organism'
        ordering=['organism_name']
        indexes = [
            models.Index(name="org_stid_idx", fields=['strain_id']),
            models.Index(name="org_stother_idx", fields=['strain_otherids']),
            models.Index(name="org_stcode_idx", fields=['strain_code']),
            models.Index(name="org_strainid_idx", fields=['strain_type']),
            models.Index(name="org_stpanel_idx", fields=['strain_panel']),
            models.Index(name="org_mta_idx", fields=['mta_status']),
            models.Index(name="org_taxid_idx", fields=['tax_id']),
            models.Index(name="org_riskgrp_idx", fields=['risk_group']),
            models.Index(name="org_pathgrp_idx", fields=['pathogen_group']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_id} ({self.strain_code})"

    #------------------------------------------------
    def str_OrganismID(self,OrganimClass,OrganismNo) -> str:
        return f"{OrganimClass}_{OrganismNo:04d}"

    #------------------------------------------------
    def find_Next_OrganismID(self,OrganimClass,OrganismClassTypes = ['GN','GP','MB','FG','MA']) -> str:
        print(f"this is find_Next_OrganismID...OrganimsClass={OrganimClass}")
        if OrganimClass in OrganismClassTypes:
            Organism_IDSq=Sequence(OrganimClass)
            Organism_nextID = next(Organism_IDSq)
            Organism_strID = self.str_OrganismID(OrganimClass,Organism_nextID)
            
            while Organism.objects.filter(organism_name=Organism_strID).exists():
                Organism_nextID = next(Organism_IDSq)
                Organism_strID = self.str_OrganismID(OrganimClass,Organism_nextID)
                print(f"name: {Organism_strID}")
            return(Organism_strID)    
        else:
            return(None)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        print("organism save model...")
        if not self.organism_id: #Object does not exists
            print("this is save from organism model...")
            self.organism_id = self.find_Next_OrganismID(str(self.organism_name.org_class.dict_value))
            if self.organism_id: 
                super().save(*args, **kwargs)
        else:
            super().save(*args, **kwargs)

    #------------------------------------------------
    def __iter__(self):
        for field in self._meta.fields:
            if field.verbose_name in ['Organism ID', 'Organism Name',  'Strain ID', 'Strain OtherID', 'Strain Code', 'Strain Notes', 'Origin']:
                yield (field.verbose_name, field.value_to_string(self))
    

#=================================================================================================
class Organism_Batch(AuditModel):
    """
    Organism/Isolate Batch Collection
    """
#-------------------------------------------------------------------------------------------------

    Choice_Dictionary = {

        'qc_status':'QC_Status',
    }

    _SEP = ':'

    orgbatch_id  = models.IntegerField(primary_key=True, blank=False, validators=[MinValueValidator(1), MaxValueValidator(10)], verbose_name = "OrgBatch ID")
    organism_id = models.ForeignKey(Organism, verbose_name = "Organism ID", on_delete=models.DO_NOTHING) 
    batch_no  = models.IntegerField(null=False, blank=False, verbose_name = "Batch No")
    batch_notes= models.CharField(blank=True, null=True, max_length=500, verbose_name = "Batch Notes")
    qc_status = models.CharField(blank=True, null=True, max_length=10, verbose_name = "QC Notes")
    qc_record = models.CharField(blank=True, null=True, max_length=10, verbose_name = "QC Records")
    supplier = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier")
    supplier_code = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier Code")
    supplier_po = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier PO")
    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date",editable=False) # Updated by Script
    stock_level = ArrayField(models.IntegerField(default=0), size=3, verbose_name = "Stock Levels", null=True, blank=False,editable=False, default=list) # Updated by Script
    biologist =models.ForeignKey(ApplicationUser, verbose_name = "Biologist", on_delete=models.DO_NOTHING) #models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")

    #------------------------------------------------
    class Meta:
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
    def find_Next_BatchNo(self,OrganismID) -> int:
        next_BatchNo = 1
        while self.objects.filter(organism_id=OrganismID, batch_no=next_BatchNo).exists():
            next_BatchNo = next_BatchNo + 1
        return(next_BatchNo)    

    #------------------------------------------------
    def str_OrgBatchID(self,OrganismID,BatchNo) -> str:
        return(f"{OrganismID}{self._SEP}{BatchNo:02d}")


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.orgbatch_id: #Object does not exists
            Next_BatchNo = self.find_Next_BatchNo(self.Organism_ID)
            if Next_BatchNo:
                self.Batch_No = Next_BatchNo
                self.orgbatch_id = self.str_OrgBatchID(self.organism_id,Next_BatchNo)
                super().save(*args, **kwargs)
        else:
            super().save(*args, **kwargs)


#=================================================================================================
class OrgBatch_Stock(AuditModel):
    """
    Stock of Organism/Isolate Batches
    
    """
#-------------------------------------------------------------------------------------------------

    Choice_Dictionary = {
        'stock_type':'Stock_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING) 
    passage_no = models.IntegerField(null=False, blank=False, verbose_name = "Passage No")
    location_freezer = models.CharField(blank=True, null=True, max_length=80, verbose_name = "Freezer")
    location_rack = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Rack")
    location_column = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Column")
    location_slot = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Slot")
    stock_type = models.CharField(blank=True, null=True, max_length=20, verbose_name = "Stock Type")
    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date")
    n_created = models.IntegerField(null=False, blank=False, verbose_name = "#Vials created")
    n_left = models.IntegerField(null=False, blank=False, verbose_name = "#Vials left")
    biologist =models.ForeignKey(ApplicationUser, verbose_name = "Biologist", on_delete=models.DO_NOTHING) # models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")

    #------------------------------------------------
    class Meta:
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
        return f"{self.orgbatch_id} {self.stock_type} {self.n_left}"

#=================================================================================================
class Organism_Culture(AuditModel):
    """
    Recommanded and optimised Growth/Culture conditions 
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'culture_type':'Culture_Type',
        'media_use':'Media_Use',
    }

    organism_id = models.ForeignKey(Organism, verbose_name = "Organism ID", on_delete=models.DO_NOTHING)
    culture_type = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Culture Type") 
    media_use = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Media Use") 
    media = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Media") 
    atmosphere = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Atmosphere") 
    temperature = models.CharField(blank=True, null=True, max_length=25, verbose_name = "Temperature") 
    labware = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Labware") 
    biologist = models.ForeignKey(ApplicationUser, verbose_name = "Biologist", on_delete=models.DO_NOTHING) #models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")
    notes = models.CharField(blank=True, null=True, max_length=512, verbose_name = "Media") 

    #------------------------------------------------
    class Meta:
        db_table = 'organism_culture'
        ordering=['organism_id','culture_type','media_use']
        indexes = [
            models.Index(name="orgcult_media_idx",fields=['media_use']),
            models.Index(name="orgcult_ctype_idx",fields=['culture_type']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organism_id} {self.media_use} {self.culture_type}"
