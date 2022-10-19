from django_rdkit import models
from django.contrib.postgres.fields import ArrayField
from app.models import Dictionaries

from sequences import Sequence
from model_utils import Choices



class Taxo(models.Model):
    """
    Based on the NCBI Taxonomy at https://www.ncbi.nlm.nih.gov/taxonomy

    with available link https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=<Tax_ID>
    NCBI has an API but with Search/Fetch procedures https://www.ncbi.nlm.nih.gov/books/NBK25501/ 

    NCBI Taxonomy specific fields (could be defined as models.TextChoice classes):
        Tax_Rank        subspecies,species,species group,genus,order,class,family,phylum,no rank,varietas
        Division        Rodents, Bacteria, Mammals, Plants and Fungi, Primates
        Division_Code   ROD, BCT, MAM, PLN, PRI

    """
 

    Name = models.CharField(primary_key=True, unique=True, max_length=150, verbose_name = "Specietest")
    
    Class_test = models.ForeignKey(Dictionaries, blank=True, null=True, verbose_name = "Class", related_name="classtest", on_delete=models.DO_NOTHING)
   
