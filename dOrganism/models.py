from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from django.contrib.postgres.fields import ArrayField
from apputil.models import AuditModel, Dictionaries
from django.db import transaction, IntegrityError

#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------


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


    def __str__(self) -> str:
        return f"{self.Organism_Name}"

    

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
    Choice_Dictionaries = {
        'Risk_Group':'Risk_Group',
        'Pathogen_Group':'Pathogen_Group',
        'Bio_Approval':'Biol_Approval',
        'Oxygen_Pref':'Oxygen_Preference',
        'MTA_Status':'License_Status',
        'Strain_Type':'Strain_Type',
        'Organism_Class':'Organism_Class',
    }


    Organism_ID = models.CharField(primary_key=True, unique=True, max_length=15, verbose_name = "Organism ID") #be blank for automatic generate a new one?
    Organism_Name= models.ForeignKey(Taxonomy, verbose_name = "Organism Name", on_delete=models.DO_NOTHING) #models do nothing?
    Organism_Desc= models.CharField(blank=True, max_length=150, verbose_name = "Organism Description", null=True)
    Strain_ID= models.CharField(blank=True, max_length=50, verbose_name = "Strain ID", null=True)
    Strain_Code= models.CharField(blank=True, max_length=15, verbose_name = "Strain Code",  null=True)
    Strain_Desc= models.CharField(blank=True, max_length=150, verbose_name = "Strain Description", null=True)
    Strain_Notes= models.CharField(blank=True, max_length=250, verbose_name = "Strain Notes", null=True)
    Strain_Tissue= models.CharField(blank=True, max_length=250, verbose_name = "Strain Tissue",  null=True)
    Strain_Type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, null=True, blank=True)
    Sequence = models.CharField(blank=True, max_length=100, verbose_name = "Sequence", null=True)
    Sequence_Link = models.CharField(blank=True, max_length=500, verbose_name = "Sequence Link", null=True)
    Tax_ID = models.IntegerField(verbose_name = "NCBI Tax ID", default=0, null=True)
    Risk_Group = models.CharField(blank=True, null=True, max_length=50)
    Pathogen_Group = models.CharField(blank=True, null=True, max_length=50)
    Import_Permit = models.CharField(blank=True, max_length=500, verbose_name = "Import Permit",  null=True)
    Bio_Approval = models.CharField(blank=True, max_length=200, verbose_name = "Biological Approval", null=True) # Dictionaries[Dictionary_ID = "Bio_Approval"]
    Special_Precaution = models.CharField(blank=True, max_length=500, verbose_name = "Special Precaution",  null=True)
    Lab_Restriction = models.CharField(blank=True, max_length=500, verbose_name = "Special Precaution",  null=True)
    MTA_Document = models.CharField(blank=True, max_length=150, verbose_name = "MTA Document", null=True)
    MTA_Status = models.CharField(blank=True, max_length=150,verbose_name = "MTA Status", null=True) # Dictionaries[Dictionary_ID = "License_Status"]
    Oxygen_Pref = models.CharField(blank=True, null=True, max_length=250)
    Atmosphere_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Atmosphere Preference", null=True)
    Nutrient_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Nutirent Preference", null=True)
    Biofilm_Pref = models.CharField(blank=True, max_length=500, verbose_name = "Biofilm Preference",null=True)

    class Meta:
        ordering=['Organism_Name']

    def __str__(self) -> str:
        return f"{self.Organism_ID} ({self.Strain_Code})"

    def save(self, *args, **kwargs):
        print('this is save function from model class')
       
        if not self.Organism_ID: #Object does not exists
            try:
                print(self.Organism_Name)
                MakeidOrganism_IDSq=Sequence(str(self.Organism_Name.Class.Dict_Value))
                MakeidOrganism_IDSq=next(MakeidOrganism_IDSq)
                self.Organism_ID=str(self.Organism_Name.Class.Dict_Value)+'_'+str(MakeidOrganism_IDSq).zfill(4)
                super().save(*args, **kwargs)
            except Exception as err:
                print(err)
                # MakeidOrganism_IDSq=Sequence('noClass')
                # MakeidOrganism_IDSq=next(MakeidOrganism_IDSq)
                # self.Organism_ID='noClass'+str(MakeidOrganism_IDSq).zfill(4)
        else:
            super().save(*args, **kwargs)

    def __iter__(self):
        for field in self._meta.fields:
            yield (field.verbose_name, field.value_to_string(self))
    


           



