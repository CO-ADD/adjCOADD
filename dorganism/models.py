from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from apputil.models import AuditModel, Dictionaries
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

    Organism_Name = models.CharField(primary_key=True, unique=True, max_length=100, verbose_name = "Specie")
    Other_Names = models.CharField(blank=True, max_length=100, verbose_name = "Other Names")
    Code = models.CharField(blank=True, max_length=15, verbose_name = "Code")
    Class = models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID")
    Parent_Tax_ID = models.IntegerField(verbose_name = "NCBI Parent Tax ID") #empty from no.207668 
    Tax_Rank = models.CharField(blank=True, max_length=50, verbose_name = "Taxonomy Rank")
    Division = models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Division", related_name='requests_Div', on_delete=models.DO_NOTHING)
    Lineage = ArrayField(models.CharField(max_length=25, null=True, blank=True),size = 25)
    
    class Meta:
        ordering=['Organism_Name']
        db_table = 'taxonomy'


    def __str__(self) -> str:
        return f"{self.Organism_Name}"

    




#-------------------------------------------------------------------------------------------------
class Organism(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionaries = {
        'Risk_Group':'Risk_Group',
        'Pathogen_Group':'Pathogen_Group',
        'Bio_Approval':'Biol_Approval',
        'Oxygen_Pref':'Oxygen_Preference',
        'MTA_Status':'License_Status',
        'Strain_Type':'Strain_Type',
        'Strain_Panel':'Strain_Panel',
        'Organism_Class':'Organism_Class',
    }


    Organism_ID = models.CharField(primary_key=True, unique=True, blank=True, max_length=15, verbose_name = "Organism ID") #be blank for automatic generate a new one?
    Organism_Name= models.ForeignKey(Taxonomy, verbose_name = "Organism Name", on_delete=models.DO_NOTHING) #models do nothing?
    #Organism_Desc= models.CharField(blank=True, max_length=150, verbose_name = "Organism Description")
    Strain_ID = models.CharField(blank=False, null=False, max_length=50, db_index=True, verbose_name = "Strain ID")
    Strain_OtherID = models.CharField(blank=True, null=True, max_length=150, verbose_name = "Strain OtherID")
    Strain_Code= models.CharField(blank=True, null=True, max_length=15, verbose_name = "Strain Code")
    #Strain_Desc= models.CharField(blank=True, max_length=150, verbose_name = "Strain Description")
    Strain_Notes= models.CharField(blank=True, null=True, max_length=250, verbose_name = "Strain Notes")
    Strain_Origin = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Origin")
    Strain_Property = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Property")
    Strain_Type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    Strain_Panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    Strain_Tissue= models.CharField(blank=True, null=True, max_length=250, verbose_name = "Strain Tissue")
    Biologist = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")
    Sequence = models.CharField(blank=True, null=True, max_length=100, verbose_name = "Sequence")
    Sequence_Link = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Sequence Link")
    Tax_ID = models.IntegerField(null=True, verbose_name = "NCBI Tax ID", default=0)
    Risk_Group = models.CharField(blank=True, null=True, max_length=50,verbose_name = "Risk Group")
    Pathogen_Group = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Pathogen",)
    Import_Permit = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Import Permit")
    Bio_Approval = models.CharField(blank=True, null=True, max_length=200, verbose_name = "Biological Approval") # Dictionaries[Dictionary_ID = "Bio_Approval"]
    Special_Precaution = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Special Precaution")
    Lab_Restriction = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Lab")
    MTA_Status = models.CharField(blank=True, null=True, max_length=150,verbose_name = "MTA Status") # Dictionaries[Dictionary_ID = "License_Status"]
    MTA_Document = models.CharField(blank=True, null=True, max_length=150, verbose_name = "MTA Document")
    Oxygen_Pref = models.CharField(blank=True, null=True, max_length=250, verbose_name = "Oxygen")
    Atmosphere_Pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Atmosphere")
    Nutrient_Pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Nutirent")
    Biofilm_Pref = models.CharField(blank=True, null=True, max_length=500, verbose_name = "Biofilm")

    class Meta:
        ordering=['Organism_Name']
        indexes = [
            models.Index(name="org_stid_idx", fields=['Strain_ID']),
            models.Index(name="org_stother_idx", fields=['Strain_OtherID']),
            models.Index(name="org_stcode_idx", fields=['Strain_Code']),
            models.Index(name="org_strainid_idx", fields=['Strain_Type']),
            models.Index(name="org_stpanel_idx", fields=['Strain_Panel']),
            models.Index(name="org_mta_idx", fields=['MTA_Status']),
            models.Index(name="org_taxid_idx", fields=['Tax_ID']),
            models.Index(name="org_riskgrp_idx", fields=['Risk_Group']),
            models.Index(name="org_pathgrp_idx", fields=['Pathogen_Group']),
        ]
        db_table = 'organism'

    def __str__(self) -> str:
        return f"{self.Organism_ID} ({self.Strain_Code})"

    def str_OrganismID(self,OrganimsClass,OrganismNo) -> str:
        return f"{OrganimsClass}_{OrganismNo:04d}"

    def find_Next_OrganismID(self,OrganimsClass,OrganismClassTypes = ['GN','GP','MB','FG','MA']) -> str:
        print(f"this is find_Next_OrganismID...OrganimsClass={OrganimsClass}")
        if OrganimsClass in OrganismClassTypes:
            Organism_IDSq=Sequence(OrganimsClass)
            Organism_nextID = next(Organism_IDSq)
            Organism_strID = self.str_OrganismID(OrganimsClass,Organism_nextID)
            
            while Organism.objects.filter(Organism_Name=Organism_strID).exists():
                Organism_nextID = next(Organism_IDSq)
                Organism_strID = self.str_OrganismID(OrganimsClass,Organism_nextID)
                print(f"name: {Organism_strID}")
            return(Organism_strID)    
        else:
            return(None)

    def save(self, *args, **kwargs):
        if not self.Organism_ID: #Object does not exists
            print("this is save from organism model...")
            self.Organism_ID = self.find_Next_OrganismID(str(self.Organism_Name.Class.Dict_Value))
            if self.Organism_ID:
                super().save(*args, **kwargs)
            
        else:
            super().save(*args, **kwargs)

    def __iter__(self):
        for field in self._meta.fields:
            yield (field.verbose_name, field.value_to_string(self))
    

#-------------------------------------------------------------------------------------------------
class Organism_Batch(AuditModel):
    """
    Organism/Isolate Batch Collection
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionaries = {
        'QC_Status':'QC_Status',
    }

    _SEP = ':'

    OrgBatch_ID  = models.IntegerField(primary_key=True, blank=False, validators=[MinValueValidator(1), MaxValueValidator(10)], verbose_name = "OrgBatch ID")
    Organism_ID = models.ForeignKey(Organism, verbose_name = "Organism ID", on_delete=models.DO_NOTHING) 
    Batch_No  = models.IntegerField(null=False, blank=False, verbose_name = "Batch No")
    Batch_Notes= models.CharField(blank=True, null=True, max_length=500, verbose_name = "Batch Notes")
    QC_Status = models.CharField(blank=True, null=True, max_length=10, verbose_name = "QC Notes")
    QC_Record = models.CharField(blank=True, null=True, max_length=10, verbose_name = "QC Records")
    Supplier = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier")
    Supplier_Code = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier Code")
    Supplier_PO = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Supplier PO")
    Stock_Date = models.DateField(null=True, blank=True, verbose_name = "Stock Date",editable=False) # Updated by Script
    Stock_Level = ArrayField(models.IntegerField(default=0), size=3, verbose_name = "Stock Levels", null=True, blank=False,editable=False, default=list) # Updated by Script
    Biologist = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")

    def find_Next_BatchNo(self,Organism_ID) -> int:
        next_BatchNo = 1
        while self.objects.filter(Organism_ID=Organism_ID, Batch_No=next_BatchNo).exists():
            next_BatchNo = next_BatchNo + 1
        return(next_BatchNo)    

    def str_OrgBatchID(self,Organism_ID,Batch_No) -> str:
        return(f"{Organism_ID}{self._SEP}{Batch_No:02d}")

    class Meta:
        ordering=['OrgBatch_ID']
        indexes = [
            models.Index(name="orgbatch_orgbatch_idx",fields=['Organism_ID','Batch_No']),
            models.Index(name="orgbatch_qc_idx",fields=['QC_Status']),
            models.Index(name="orgbatch_supp_idx",fields=['Supplier']),
            models.Index(name="orgbatch_sdate_idx",fields=['Stock_Date']),
            models.Index(name="orgbatch_slevel_idx",fields=['Stock_Level']),
        ]
        db_table = 'organism_batch'

    def save(self, *args, **kwargs):
        if not self.OrgBatch_ID: #Object does not exists
            Next_BatchNo = self.find_Next_BatchNo(self.Organism_ID)
            if Next_BatchNo:
                self.Batch_No = Next_BatchNo
                self.OrgBatch_ID = self.str_OrgBatchID(self.Organism_ID,Next_BatchNo)
                super().save(*args, **kwargs)
        else:
            super().save(*args, **kwargs)

    def __str__(self) -> str:
        return f"{self.OrgBatch_ID}"

#-------------------------------------------------------------------------------------------------
class OrgBatch_Stock(AuditModel):
    """
    Stock of Organism/Isolate Batches
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionaries = {
        'Stock_Type':'Stock_Type',
    }

    OrgBatch_ID = models.ForeignKey(Organism_Batch, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING) 
    Passage_No = models.IntegerField(null=False, blank=False, verbose_name = "Passage No")
    Location_Freezer = models.CharField(blank=True, null=True, max_length=80, verbose_name = "Freezer")
    Location_Rack = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Rack")
    Location_Column = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Column")
    Location_Slot = models.CharField(blank=True, null=True, max_length=10, verbose_name = "Slot")
    Stock_Type = models.CharField(blank=True, null=True, max_length=20, verbose_name = "Stock Type")
    Stock_Date = models.DateField(null=True, blank=True, verbose_name = "Stock Date")
    N_Created = models.IntegerField(null=False, blank=False, verbose_name = "#Vials created")
    N_Left = models.IntegerField(null=False, blank=False, verbose_name = "#Vials left")
    Biologist = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")

    class Meta:
        ordering=['OrgBatch_ID','Stock_Type']
        indexes = [
            models.Index(name="orgbstock_stype_idx",fields=['Stock_Type']),
            models.Index(name="orgbstock_freezer_idx",fields=['Location_Freezer']),
            models.Index(name="orgbstock_stdate_idx",fields=['Stock_Date']),
            models.Index(name="orgbstock_nleft_idx",fields=['N_Left']),
        ]
        db_table = 'orgbatch_stock'
    def __str__(self) -> str:
        return f"{self.OrgBatch_ID} {self.Stock_Type} {self.N_Left}"

#-------------------------------------------------------------------------------------------------
class Organism_Culture(AuditModel):
    """
    Recommanded and optimised Growth/Culture conditions 
    
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionaries = {
        'Culture_Type':'Culture_Type',
        'Media_Use':'Media_Use',
    }

    Organism_ID = models.ForeignKey(Organism, verbose_name = "Organism ID", on_delete=models.DO_NOTHING)
    Culture_Type = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Culture Type") 
    Media_Use = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Media Use") 
    Media = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Media") 
    Atmosphere = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Atmosphere") 
    Temperature = models.CharField(blank=True, null=True, max_length=25, verbose_name = "Temperature") 
    Labware = models.CharField(blank=True, null=True, max_length=120, verbose_name = "Labware") 
    Biologist = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Biologist")
    Notes = models.CharField(blank=True, null=True, max_length=512, verbose_name = "Media") 

    class Meta:
        ordering=['Organism_ID','Culture_Type','Media_Use']
        indexes = [
            models.Index(name="orgcult_media_idx",fields=['Media_Use']),
            models.Index(name="orgcult_ctype_idx",fields=['Culture_Type']),
        ]
        db_table = 'organism_culture'
    def __str__(self) -> str:
        return f"{self.Organism_ID} {self.Media_Use} {self.Culture_Type}"