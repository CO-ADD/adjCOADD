from django.db import models
import re
from model_utils import Choices
from sequences import Sequence
#from django_rdkit import models
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from django.core.validators import RegexValidator
from dorganism.models import Taxonomy 

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Cell Application Model
#-------------------------------------------------------------------------------------------------

class Cell(AuditModel):

    cell_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Cell ID") 
    cell_names= models.CharField(max_length=200, blank=True, verbose_name = "Cell Names") 
    cell_line= models.CharField(max_length=200, blank=True, verbose_name = "Cell Line") 
    cell_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Cell Notes")
    cell_code= models.CharField(max_length=30, blank=True, verbose_name = "Cell Code")
    cell_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    cell_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    cell_identification = models.CharField(max_length=512, blank=True, verbose_name = "Cell Identification")
    cell_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Cell")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
        db_column="organism_name", related_name="%(class)s_organism_name")

    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_mta")
    mta_notes = models.CharField(max_length=150, blank=True, verbose_name = "MTA Notes")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    mta_notes = models.CharField(max_length=512, blank=True, verbose_name = "MTA Notes")

    collect_tissue = models.CharField(max_length=120, blank=True, verbose_name = "From Tissue/Organ")
    patient_diagnosis = models.CharField(max_length=120, blank=True, verbose_name = "Patient Diagnosis")
    patient = models.CharField(max_length=20, blank=True, verbose_name = "Patient Info")

    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    assoc_documents = models.ManyToManyField(Document,verbose_name = "Douments", blank=True,
        db_table = "org_doc", related_name="%(class)s_document")
    
    #tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")


#------------------------------------------------


    class Meta:
        app_label = 'dcell'
        db_table = 'cell'
        #ordering=['organism_id']
        indexes = [
            models.Index(name="cell_line_idx", fields=['cell_line']),
            # models.Index(name="org_stcode_idx", fields=['strain_code']),
            # models.Index(name="org_strainid_idx", fields=['strain_type']),
            # models.Index(name="org_stpanel_idx", fields=['strain_panel']),
            # models.Index(name="org_source_idx", fields=['source']),
            # models.Index(name="org_taxid_idx", fields=['tax_id']),
            # models.Index(name="org_riskgrp_idx", fields=['risk_group']),
            # models.Index(name="org_pathgrp_idx", fields=['pathogen_group']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.cell_id} ({self.cell_line})"

    #------------------------------------------------
    @classmethod
    def str_CellID(cls,CellNo) -> str:
    #
    # Input:    OrganismClass GN, GP,...
    #           OrganismNo 
    # Output:   Oragnism_ID as string like GN_0001 
    #
        return(f"CL{ORGANSIM_SEP}{CellNo:04d}")

    #------------------------------------------------
    @classmethod
    def exists(cls,CellID,verbose=0):
    # Returns if an instance exists by organism_id
        return cls.objects.filter(cell_id=CellID).exists()

    #------------------------------------------------
    @classmethod
    def get(cls,CellID,verbose=0):
    # Returns an instance by organism_id
        try:
            retInstance = cls.objects.get(cell_id=CellID)
        except:
            if verbose:
                print(f"[CellID Not Found] {CellID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def find_Next_OrganismID(cls) -> str:
        Cell_IDSq = Sequence("Cell")
        Cell_nextID = next(Cell_IDSq)
        Cell_strID = cls.str_CellID(Cell_nextID)
        while cls.exists(Cell_strID):
            Cell_nextID = next(Cell_IDSq)
            Cell_strID = cls.str_CellID(Cell_nextID)
        return(Cell_strID)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.cell_id: 
            self.cell_id = self.find_Next_CellID()
            if self.cell_id: 
                super(Cell, self).save(*args, **kwargs)
        else:
            super(Cell, self).save(*args, **kwargs) 
