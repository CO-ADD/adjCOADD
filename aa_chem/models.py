from django_rdkit import models
from django.contrib.postgres.fields import ArrayField
from app.models import User
# Create your models here.

class Drugbank(models.Model):
    
    status = models.IntegerField(blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True,blank=True, null=True)
    modified_at = models.DateTimeField(auto_now=True, blank=True, null=True)
    drug_id = models.CharField(unique=True, max_length=15, blank=True, null=True)
    drug_name = models.CharField(max_length=50, blank=True, null=True)
    drug_code = models.CharField(max_length=10, blank=True, null=True)
    drug_synonyms = models.CharField(max_length=1024, blank=True, null=True)
    synonyms = models.TextField(blank=True, null=True)  # This field type is a guess.
    cas = models.CharField(max_length=15, blank=True, null=True)
    unii = models.CharField(max_length=15, blank=True, null=True)
    access_ids = models.CharField(max_length=250, blank=True, null=True)
    db_access_ids = models.TextField(blank=True, null=True)  # This field type is a guess.
    drug_mol = models.MolField(blank=True, null=True)  # This field type is a guess.
    drug_smiles = models.CharField(max_length=2048, blank=True, null=True)
    
    def save(self, *args, **kwargs):
        if not self.drug_id:
            self.drug_id='us'+str(Drugbank.objects.count()+1)  #drug_id is usXXXXX meaning...
            
        super().save(*args, **kwargs)


    class Meta:
        managed = True
        db_table = 'drugbank'



#
# Organism/Strain Models
#

from typing import Sequence
# from django.db import models
from model_utils import Choices

class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
    astatus = models.IntegerField(verbose_name = "Status", default = 0, editable=False) #, index = True
    acreated_at = models.DateTimeField(auto_now_add=True, null=True, verbose_name = "Created at",editable=False)
    aupdated_at = models.DateTimeField(auto_now=True, null=True, blank=True, verbose_name = "Updated at",editable=False)
    adeleted_at = models.DateTimeField(blank=True, null=True, verbose_name = "Deleted at",editable=False)
    acreated_by = models.ForeignKey(User,blank=True, null=True,editable=False, verbose_name = "Created by", related_name='%(class)s_requests_created', on_delete=models.CASCADE) #
    aupdated_by = models.ForeignKey(User, blank=True, null=True,editable=False, verbose_name = "Updated by", related_name='%(class)s_requests_updated',on_delete=models.CASCADE) #
    adeleted_by = models.ForeignKey(User, blank=True, null=True,editable=False, verbose_name = "Deleted by", related_name='%(class)s_requests_deleted',on_delete=models.CASCADE) #

    class Meta:
        abstract = True
        indexes = [
            models.Index(fields=['astatus']),
            
        ]
    
    def delete(self,**kwargs):
        self.astatus = -9
        self.adeleted_at = timezone.now()
        self.adeleted_by = kwargs.get("user")
        self.save(updated_fields = ['adeleted_at','adeleted_by','astatus'])
    
    def save(self, *args, **kwargs):
        user = kwargs.get("user")
        if self.pk: #Object already exists
            self.aupdated_by = user
        else:
            self.acreated_by = user
        super(AuditModel,self).save(*args, **kwargs)


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

    """

    Divisions = Choices(
        ('ROD', ('Rodents')), 
        ('BCT', ('Bacteria')),
        ('MAM', ('Mammals')),
        ('PLN', ('Plants and Fungi')),
        ('PRI', ('Primates')),
    )

    Classes = Choices(
        ('GN', ('Gram-Negative')), 
        ('GP', ('Gram-Positive')),
        ('MB', ('Mycobacteria')),
        ('FG', ('Fungi')),
        ('PL', ('Plant')),
        ('MA', ('Mammalian')),
    )

    Organism_Name = models.CharField(unique=True, max_length=150, verbose_name = "Specie")
    Other_Names = models.CharField(blank=True, max_length=250, verbose_name = "Other Names")
    Code = models.CharField(blank=True, max_length=12, verbose_name = "Code")
    Class = models.CharField(blank=True, max_length=50, verbose_name = "Class",choices=Classes)
    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID")
    Parent_Tax_ID = models.IntegerField(verbose_name = "NCBI Parent Tax ID", default=468123) #empty from no.207668 
    Tax_Rank = models.CharField(blank=True, max_length=50, verbose_name = "Taxonomy Rank")
    Division = models.CharField(blank=True, max_length=25, verbose_name = "Division", choices=Divisions)
    Division_CODE = models.CharField(blank=True, max_length=25, verbose_name = "Division_code",default='PLN')
    
    #Lineage = models.CharField(blank=True, max_length=1024, verbose_name = "Lineage")
    

    Lineage = ArrayField( models.CharField(max_length=25, blank=True),size = 15) 
   
    
    class Meta:
        db_table = 'strain_taxonomy'

    def __str__(self) -> str:
        return f"{self.Organism_Name}"

class Genes(AuditModel):
    """
    List of different bacterial genes  of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    """

#-------------------------------------------------------------------------------------------------
"""
Sequence numbers https://pypi.org/project/django-sequences/
"""
from sequences import Sequence

GN_Sequence=Sequence("Gram-Negative")
GP_Sequence=Sequence("Gram-Positive")
MB_Sequence=Sequence("Mycobacteria")
FG_Sequence=Sequence("Fungi")
MA_Sequence=Sequence("Mammalian")

class Organisms(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """

    RiskGroups = Choices(
        ('RG1', ('Risk Group 1')),
        ('RG2', ('Risk Group 2')),
    )

    PathogenGroups = Choices(
        ('PAT', ('Pathogen')),
        ('NOP', ('Non-Pathogen')),
    )

    Organism_ID = models.CharField(unique=True, blank=True, max_length=15, verbose_name = "OrgID") #must be blank for automatic generate a new one
    Organism_Name_set= models.ForeignKey(Taxonomy, verbose_name = "Organism Name Set",on_delete=models.DO_NOTHING) #models do nothing?
    Organism_Name=models.CharField(unique=True,blank=True, max_length=312, verbose_name = "Organism Name", editable=False)
    Organism_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Organism Description", default="NotFoundValue")
    # Organism_Class= models.CharField(blank=True, max_length=50, verbose_name = "Organism Class") ... in Taxonomy.Organism_Class

    Strain_ID= models.CharField(blank=True, max_length=250, verbose_name = "Strain ID", default="NotFoundValue")
    Strain_Code= models.CharField(blank=True, max_length=500, verbose_name = "Strain Code", default="NotFoundValue")
    Strain_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Strain Description", default="NotFoundValue")
    Strain_Notes= models.CharField(blank=True, max_length=512, verbose_name = "Strain Notes", default="NotFoundValue")
    Strain_Tissue= models.CharField(blank=True, max_length=220, verbose_name = "Strain Tissue", default="NotFoundValue")
    Strain_Type= models.CharField(blank=True, max_length=250, verbose_name = "Strain Types", default="NotFoundValue")

    Sequence = models.CharField(blank=True, max_length=512, verbose_name = "Sequence", default="NotFoundValue")
    Sequence_Link = models.CharField(blank=True, max_length=1000, verbose_name = "Sequence Link", default="NotFoundValue")
    Geno_Type = models.CharField(blank=True, max_length=512, verbose_name = "GenoType", default="NotFoundValue")

    Screen_Type = ArrayField(
        models.CharField(max_length=1000, blank=True, null=True), size=20, null=True, blank=True
    )

    """
    # Optional for Strain_Type as Array based on multi-selectable
    #  using LineageArray = Strain_Type.split('; ')
    Strain_Type = models.ArrayField( models.CharField(max_length=10, blank=True),size = 15) 
    """

    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID", default=0)

    Risk_Group = models.CharField(blank=True, max_length=500, verbose_name = "Risk Group",choices=RiskGroups, default="NotFoundValue")
    Pathogen = models.CharField(blank=True, max_length=500, verbose_name = "Risk Group",choices=PathogenGroups, default="NotFoundValue")

    Import_Permit = models.CharField(blank=True,max_length=500, verbose_name = "Import Permit", default="NotFoundValue")
    Biol_Approval = models.CharField(blank=True, max_length=500, verbose_name = "Biological Approval", default="NotFoundValue")
    Special_Precaution = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution", default="NotFoundValue")
    Lab_Restriction = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution", default="NotFoundValue")
    MTA_Document = models.CharField(blank=True, max_length=500, verbose_name = "MTA Document", default="NotFoundValue")
    MTA_Status = models.CharField(blank=True, max_length=500, verbose_name = "MTA Status", default="NotFoundValue")

    Oxygen_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Oxygen Preference")
    Atmosphere_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Atmosphere Preference")
    Nutrient_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Nutirent Preference")
    Biofilm_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Biofilm Preference")

    class Meta:
        db_table = 'strain_organisms'

    def __str__(self) -> str:
        return f"{self.Organism_ID} ({self.Strain_Code})"

    def save(self, *args, **kwargs):
        
        if not self.Organism_ID: #Object does not exists
            org_id_pre=Sequence(self.Organism_Name_set.Class)

            num=next(org_id_pre)
            self.Organism_ID=self.Organism_Name_set.Class+'_'+'000'+str(num)
            # print(self.Organism_ID)
            # self.acreated_by = user
        if not self.Organism_Name:
            self.Organism_Name=self.Organism_Name_set.Organism_Name
        super().save(*args, **kwargs)

