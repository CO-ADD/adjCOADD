from django_rdkit import models
from django.contrib.postgres.fields import ArrayField

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
from django.db import models
from model_utils import Choices

class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
    astatus = models.IntegerField(verbose_name = "Status", default = 0, index = True, editable=False)
    acreated_at = models.DateTimeField(auto_now_add=True, verbose_name = "Created at",editable=False)
    aupdated_at = models.DateTimeField(auto_now=True, blank=True, verbose_name = "Updated at",editable=False)
    adeleted_at = models.DateTimeField(blank=True, verbose_name = "Deleted at",editable=False)
    acreated_by = models.ForeignKey(User, verbose_name = "Created by",editable=False)
    aupdated_by = models.ForeignKey(User, blank=True, verbose_name = "Updated by",editable=False)
    adeleted_by = models.ForeignKey(User, blank=True, verbose_name = "Deleted by",editable=False)

    class Meta:
        abstract = True
    
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
        ('ROD',_('Rodents')), 
        ('BCT',_('Bacteria')),
        ('MAM',_('Mammals')),
        ('PLN',_('Plants and Fungi')),
        ('PRI',_('Primates')),
    )

    Classes = Choices(
        ('GN',_('Gram-Negative')), 
        ('GP',_('Gram-Positive')),
        ('MB',_('Mycobacteria')),
        ('FG',_('Fungi')),
        ('PL',_('Plant')),
        ('MA',_('Mammalian')),
    )

    Organism_Name = models.CharField(unique=True, max_length=150, verbose_name = "Specie")
    Other_Names = models.CharField(blank=True, max_length=250, verbose_name = "Other Names")
    Code = models.CharField(blank=True, max_length=12, verbose_name = "Code")
    Class = models.CharField(blank=True, max_length=50, verbose_name = "Class",choices=Classes)
    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID")
    Parent_Tax_ID = models.IntegerField(verbose_name = "NCBI Parent Tax ID")
    Tax_Rank = models.CharField(blank=True, max_length=50, verbose_name = "Taxonomy Rank")
    Division = models.CharField(blank=True, max_length=25, verbose_name = "Division",choices=Divisions)
    Lineage = models.CharField(blank=True, max_length=1024, verbose_name = "Lineage")
    """
    # Optional for Linage as Array 
    #  using LineageArray = Lineage.split('; ')

    Lineage = models.ArrayField( models.CharField(max_length=25, blank=True),size = 15) 
    """
    
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

from sequences import Sequence

GN_Sequence=Sequence("Gram-Negative")
GP_Sequence=Sequence("Gram-Positive")
MB_Sequence=Sequence("Mycobacteria")
FG_Sequence=Sequence("Fungi")
MA_Sequence=Sequence("Mammalian")
"""
class Organisms(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """

    RiskGroups = Choices(
        ('RG1',_('Risk Group 1')),
        ('RG2',_('Risk Group 2')),
    )

    PathogenGroups = Choices(
        ('PAT',_('Pathogen')),
        ('NOP',_('Non-Pathogen')),
    )

    Organism_ID = models.CharField(unique=True, max_length=15, verbose_name = "OrgID")
    Organism_Name = models.ForeignKey(Taxonomy, blank=True, verbose_name = "Organism Name",on_delete=models.DO_NOTHING)
    Organism_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Organism Description")
    # Organism_Class= models.CharField(blank=True, max_length=50, verbose_name = "Organism Class") ... in Taxonomy.Organism_Class

    Strain_ID= models.CharField(blank=True, max_length=150, verbose_name = "Strain ID")
    Strain_Code= models.CharField(blank=True, max_length=50, verbose_name = "Strain Code")
    Strain_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Strain Description")
    Strain_Notes= models.CharField(blank=True, max_length=512, verbose_name = "Strain Notes")
    Strain_Tissue= models.CharField(blank=True, max_length=120, verbose_name = "Strain Tissue")
    Strain_Type= models.CharField(blank=True, max_length=150, verbose_name = "Strain Types")

    Sequence = models.CharField(blank=True, max_length=512, verbose_name = "Sequence")
    Sequence_Link = models.CharField(blank=True, max_length=120, verbose_name = "Sequence Link")
    Geno_Type = models.CharField(blank=True, max_length=512, verbose_name = "GenoType")

    Screen_Type = 

    """
    # Optional for Strain_Type as Array based on multi-selectable
    #  using LineageArray = Strain_Type.split('; ')
    Strain_Type = models.ArrayField( models.CharField(max_length=10, blank=True),size = 15) 
    """

    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID")

    Risk_Group = models.CharField(blank=True, max_length=3, verbose_name = "Risk Group",choices=RiskGroups)
    Pathogen = models.CharField(blank=True, max_length=3, verbose_name = "Risk Group",choices=PathogenGroups)

    Import_Permit = models.CharField(blank=True, max_length=150, verbose_name = "Import Permit")
    Biol_Approval = models.CharField(blank=True, max_length=10, verbose_name = "Biological Approval")
    Special_Precaution = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution")
    Lab_Restriction = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution")
    MTA_Document = models.CharField(blank=True, max_length=50, verbose_name = "MTA Document")
    MTA_Status = models.CharField(blank=True, max_length=50, verbose_name = "MTA Status")

    Oxygen_Pref = models.CharField(blank=True, max_length=50, verbose_name = "Oxygen Preference")
    Atmosphere_Pref = models.CharField(blank=True, max_length=50, verbose_name = "Atmosphere Preference")
    Nutrient_Pref = models.CharField(blank=True, max_length=50, verbose_name = "Nutirent Preference")
    Biofilm_Pref = models.CharField(blank=True, max_length=50, verbose_name = "Biofilm Preference")

    class Meta:
        db_table = 'strain_organisms'

    def __str__(self) -> str:
        return f"{self.OrgID} ({self.Strain_Code})"

    def save(self, *args, **kwargs):
        if ~self.pk: #Object does not exists
            if self.
             self.acreated_by = user
        self.save(*args, **kwargs)

