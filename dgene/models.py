from django.db import models

from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary
from dorganism.models import Organism, Organism_Batch
from dscreen.models import Screen_Run

#-------------------------------------------------------------------------------------------------
# Gene and Identificationrelated Application Model
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Gene(AuditModel):
    """
    List of Genes
    """
#=================================================================================================
    HEADER_FIELDS = {
        "gene_id":"Gene ID",
        "gene_name":"Gene Name",
        "gene_othernames":"Other Name",
        "gene_codes":"Gene Code",
        "gene_type":"Gene Type",
        "gene_class":"Gene Class",
        "gene_note":"Gene Note",
    }

    Choice_Dictionary = {
        'gene_type':'Gene_Type',
        'max_phase':'Max_Phase',
    }

#=================================================================================================
class ID_Pub(AuditModel):
    """
     Identification from Public or Collaborative sources    
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "organism_id":"Organism ID",
        "id_type":"ID Method",
        "id_organism":"Organism",
        "id_probability": "Probability",
        "id_confidence": "Confidence",
        "id_source": "Source",
        "id_date":"Date",
    }
    Choice_Dictionary = {
        'id_type':'ID_Type',
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id") 

    id_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "ID Method", on_delete=models.DO_NOTHING,
         db_column="id_type", related_name="%(class)s_idtype")
    id_organism = models.CharField(max_length=120,  blank=True, verbose_name = "ID Organism")
    id_probability = models.CharField(max_length=120, blank=True,  verbose_name = "ID Probability")
    id_confidence = models.CharField(max_length=120, blank=True,  verbose_name = "ID Confidence")
    id_source = models.CharField(max_length=20,  blank=True, verbose_name = "ID Source")
    id_date = models.DateField(blank=True, verbose_name = "ID Date")

#=================================================================================================
class ID_COADD(AuditModel):
    """
     Identification from Internal experiments    
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "orgbatch_id":"OrgBatch ID",
        "id_type":"ID Method",
        "id_organism":"Organism",
        "id_probability": "Probability",
        "id_confidence": "Confidence",
        "id_source": "Source",
        "id_date":"Date",
    }
    Choice_Dictionary = {
        'id_type':'ID_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 

    id_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "ID Method", on_delete=models.DO_NOTHING,
         db_column="id_type", related_name="%(class)s_idtype")
    id_organism = models.CharField(max_length=120,  blank=True, verbose_name = "ID Organism")
    id_probability = models.CharField(max_length=120, blank=True,  verbose_name = "ID Probability")
    id_confidence = models.CharField(max_length=120, blank=True,  verbose_name = "ID Confidence")
    id_date = models.DateField(blank=True, verbose_name = "ID Date")
    id_notes = models.CharField(max_length=120, blank=True,  verbose_name = "ID Notes")
    run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID")
