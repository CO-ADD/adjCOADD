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