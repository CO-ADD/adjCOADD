from django_rdkit import models
from django.contrib.postgres.fields import ArrayField
from app.models import AuditModel, Dictionaries, ChoiceArrayField#, Choice_Dictionaries
from sequences import Sequence
from model_utils import Choices

# from multiselectfield import MultiSelectField

# MY_CHOICES = (('item_key1', 'Item title 1.1'),
#               ('item_key2', 'Item title 1.2'),
#               ('item_key3', 'Item title 1.3'),
#               ('item_key4', 'Item title 1.4'),
#               ('item_key5', 'Item title 1.5'))

# Create your models here.

#-------------------------------------------------------------------------------------------------
class Choices_multi(models.Model):
    Unit1= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    Unit2= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    Unit3= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
#-------------------------------------------------------------------------------------------------
# Unit1='a'
# Unit2='b'
# Unit3='c'

# Choice_Dictionaries = {
#         (Choices_multi.Unit1,'Unit1'),
#         (Choices_multi.Unit1,'Unit2'),
#         (Choices_multi.Unit1,'Unit3'),
#     }

#-------------------------------------------------------------------------------------------------
class Drugbank(models.Model):
#-------------------------------------------------------------------------------------------------    
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
    # Unit1= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    # Unit2= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    # Unit3= models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    Unit1='PF'
    Unit2='F'
    Unit3='P'
    
    Choice_Dictionaries = {
        (Unit1,'Plants and Fungi'),
        (Unit2 ,'Fungi'),
        (Unit3 ,'Plants'),
    }

    Organism_Name = models.CharField(primary_key=True, unique=True, max_length=150, verbose_name = "Specie")
    Other_Names = models.CharField(blank=True, max_length=250, verbose_name = "Other Names")
    Code = models.CharField(blank=True, max_length=12, verbose_name = "Code")
    Class = models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="dictClass+", on_delete=models.DO_NOTHING)
    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID")
    Parent_Tax_ID = models.IntegerField(verbose_name = "NCBI Parent Tax ID") #empty from no.207668 
    Tax_Rank = models.CharField(blank=True, max_length=50, verbose_name = "Taxonomy Rank")
    # Division= ChoiceArrayField(models.CharField(max_length=150,choices=Choice_Dictionaries), default=list)
    Division = models.ForeignKey(Dictionaries, blank=True, verbose_name = "Division", related_name='%(class)s_requests_Div',on_delete=models.DO_NOTHING)
    Lineage = ArrayField(models.CharField(max_length=25, null=True, blank=True),size = 15)

    def __str__(self) -> str:
        return f"{self.Organism_Name}"

    def save(self, *args, **kwargs):

        super().save(*args, **kwargs)

class Genes(models.Model):
    pass
    """
    List of different bacterial genes  of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    """

#-------------------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------------------
class Organisms(AuditModel):
    """
    Main class of Organisms/Bacterias/Fungi/Cells in Isolate Collection
    
    """
#-------------------------------------------------------------------------------------------------
    # Choice_Dictionaries = {
    #     'Risk_Group':'Risk_Group',
    #     'Pathogen_Group':'Pathogen_Group',
    #     'Bio_Approval':'Bio_Approval',
    #     'Oxygen_Pref':'Oxygen_Preference',
    #     'MTA_Status':'License_Status',
    #     'Strain_Type':'Strain_Type',
    # }

    Unit1='RM'
    Unit2='CI'
    Unit3='RX'
    
    Choice_Dictionaries = {
        (Unit1,'Resistant MDR'),
        (Unit2,'Clinical Isolate'),
        (Unit3,'Resistant XDR'),
    }

    Organism_ID = models.CharField(primary_key=True, unique=True, blank=True, max_length=100, verbose_name = "Organism ID") #be blank for automatic generate a new one?
    Organism_Name= models.ForeignKey(Taxonomy, null=True, blank=True, verbose_name = "Organism Name", on_delete=models.DO_NOTHING ) #models do nothing?
    Organism_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Organism Description", default="--", null=True)
    Strain_ID= models.CharField(blank=True, max_length=250, verbose_name = "Strain ID", null=True)
    Strain_Code= models.CharField(blank=True, max_length=500, verbose_name = "Strain Code", default="--", null=True)
    Strain_Desc= models.CharField(blank=True, max_length=512, verbose_name = "Strain Description", default="--",null=True)
    Strain_Notes= models.CharField(blank=True, max_length=512, verbose_name = "Strain Notes", default="--",null=True)
    Strain_Tissue= models.CharField(blank=True, max_length=220, verbose_name = "Strain Tissue", default="--", null=True)
    Strain_Type= models.CharField(max_length=150, null=True, blank=True)
    # Strain_Type=ChoiceArrayField(models.CharField(max_length=150,choices=[]), default=list)
    # Strain_Type=models.ManyToManyField(Dictionaries)
    Sequence = models.CharField(blank=True, max_length=512, verbose_name = "Sequence", default="--",null=True)
    Sequence_Link = models.CharField(blank=True, max_length=1000, verbose_name = "Sequence Link", default="--",null=True)

    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID", default=0, null=True)

    Risk_Group = models.ForeignKey(Dictionaries,blank=True, verbose_name = "Risk Group",related_name='%(class)s_requests_RG', on_delete=models.DO_NOTHING, null=True)   # Dictionaries[Dictionary_ID = "Risk_Group"]
    Pathogen = models.ForeignKey(Dictionaries, blank=True, verbose_name = "Pathogen Group",related_name='%(class)s_requests_PT', on_delete=models.DO_NOTHING, null=True) # Dictionaries[Dictionary_ID = "Pathogen_Group"]

    Import_Permit = models.CharField(blank=True, max_length=500, verbose_name = "Import Permit", default="--", null=True)
    Biol_Approval = models.ForeignKey(Dictionaries, blank=True, verbose_name = "Biological Approval",related_name='%(class)s_requests_BA', on_delete=models.DO_NOTHING,null=True) # Dictionaries[Dictionary_ID = "Bio_Approval"]
    Special_Precaution = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution", default="--", null=True)
    Lab_Restriction = models.CharField(blank=True, max_length=512, verbose_name = "Special Precaution", default="--", null=True)
    MTA_Document = models.CharField(blank=True, max_length=500, verbose_name = "MTA Document", default="--",null=True)
    MTA_Status = models.ForeignKey(Dictionaries,blank=True, verbose_name = "MTA Status",related_name='%(class)s_requests_MTA', on_delete=models.DO_NOTHING, null=True) # Dictionaries[Dictionary_ID = "License_Status"]

    Oxygen_Pref = models.ForeignKey(Dictionaries,blank=True, verbose_name = "Oxygen Preference",related_name='%(class)s_requests_OP', on_delete=models.DO_NOTHING, null=True) # Dictionaries[Dictionary_ID = "Oxygen_Preference"]
    Atmosphere_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Atmosphere Preference", null=True)
    Nutrient_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Nutirent Preference", null=True)
    Biofilm_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Biofilm Preference",null=True)

    def __str__(self) -> str:
        return f"{self.Organism_ID} ({self.Strain_Code})"

    def save(self, *args, **kwargs):
       
        if not self.Organism_ID: #Object does not exists
            num=Sequence(str(self.Organism_Name.Class.Dict_Value))
            try:
                
                num=next(num)
                self.Organism_ID=str(self.Organism_Name.Class.Dict_Value)+'_'+str(num).zfill(4)
            except Exception as err:
                print(err)
                self.Organism_ID='__'+'_'.zfill(4)
            super().save(*args, **kwargs)
           
class Mytest(models.Model):
    
    tag=models.CharField(max_length=155, default='default') 
    date=models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.tag


