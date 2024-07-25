import re

#from django_rdkit import models
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from dchem.models import Chem_Structure
from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Screening Application Model
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
class Screen_Run(AuditModel):
    """
    List of Screening runs
    """
#-------------------------------------------------------------------------------------------------
    HEADER_FIELDS = {
        "run_id":"Run ID",
        "run_type":"Run Type",
        "assay_note":"Assay",
        "run_status":"Status",
        "run_project":"Project",
        "run_name":"Name",
        "run_date":"Run Date",
        "run_conditions":"Conditions",
        "run_issues":"Issues",
    }

    Choice_Dictionary = {
        'run_type':'Run_Type',
        'run_status':'Run_Status',
    }

    run_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Run ID")
    run_name = models.CharField(max_length=500, verbose_name = "Run Name")
    run_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Type", on_delete=models.DO_NOTHING,
        db_column="run_type", related_name="%(class)s_RunType+")
    assay_note = models.CharField(max_length=250, blank=True, verbose_name = "Assay Note")
    run_conditions = models.CharField(max_length=250, blank=True, verbose_name = "Run Conditions")
    run_issues = models.CharField(max_length=250, blank=True, verbose_name = "Run Issues")
    run_date = models.DateField(null=True, blank=True, verbose_name = "Run Date")
    run_project = models.CharField(max_length=50, verbose_name = "Project")
    run_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Status", on_delete=models.DO_NOTHING,
        db_column="run_status", related_name="%(class)s_run_status+")

    #------------------------------------------------
    class Meta:
        app_label = 'dscreen'
        db_table = 'screen_run'
        ordering=['run_type','run_id']
        indexes = [
            models.Index(name="run_type_idx", fields=['run_type']),
        ]

    # #------------------------------------------------
    # def __str__(self) -> str:
    #     return f"{self.run_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.run_id} [{self.run_type}]"

    #------------------------------------------------
    @classmethod
    def get(cls,RunID,verbose=0):
        try:
            retInstance = cls.objects.get(run_id=RunID)
        except:
            if verbose:
                print(f"[Run_ID Not Found] {RunID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,RunID,verbose=0):
        return cls.objects.filter(run_id=RunID.strip()).exists()

#-------------------------------------------------------------------------------------------------
# Samples 
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
# class Screen_Sample(AuditModel):
#     """
#     List of Compounds srceened
#     """
# #-------------------------------------------------------------------------------------------------

#     ID_SEQUENCE = 'CastDB_Sample'
#     ID_PREFIX = 'S'
#     ID_PAD = 9

#     sample_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Sample ID")
#     sample_code = models.CharField(max_length=15, blank=True, verbose_name = "Sample Code")
#     sample_name = models.CharField(max_length=50, blank=True, verbose_name = "Sample Name")
#     sample_desc = models.CharField(max_length=512, blank=True, verbose_name = "Sample Name")
#     sample_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Sample Type", on_delete=models.DO_NOTHING,
#         db_column="sample_type", related_name="%(class)s_sampletype")
#     old_compound_id = models.CharField(max_length=15, unique=True, verbose_name = "Compound ID")
#     previous_ids = models.CharField(max_length=100, blank=True, verbose_name = "Previous IDs")
#     parent_structure_ids = ArrayField(models.CharField(max_length=15, null=True, blank=True), size=4, verbose_name = "Panel", 
#                                       null=True, blank=True)
    
#     project_id = models.CharField(max_length=15, blank=True, verbose_name = "Project ID")

#     reg_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Smiles")
#     reg_mw = models.FloatField(default=0, blank=True, verbose_name = "Reg MW")
#     reg_mf = models.CharField(max_length=100, blank=True, verbose_name = "Reg MF")
#     reg_structure = models.CharField(max_length=2048, blank=True, verbose_name = "Reg Structure")
#     reg_amount = models.FloatField(default=0, blank=True, verbose_name = "Reg Amount")
#     reg_amount_unit = models.CharField(max_length=100, blank=True, verbose_name = "Reg Amount Unit")
#     reg_volume = models.FloatField(default=0, blank=True, verbose_name = "Reg Volume")
#     reg_volume_unit = models.CharField(max_length=100, blank=True, verbose_name = "Reg Volume Unit")
#     reg_conc = models.FloatField(default=0, blank=True, verbose_name = "Reg Conc")
#     reg_conc_unit = models.CharField(max_length=100, blank=True, verbose_name = "Reg Conc Unit")
#     reg_solvent = models.FloatField(default=0, blank=True, verbose_name = "Reg Solvent")
    
#     prep_date = models.DateField(null=False, editable=False, verbose_name="Prepared")
#     stock_amount = models.FloatField(default=0, blank=True, verbose_name = "Stock Amount")
#     stock_amount_unit = models.CharField(max_length=100, blank=True, verbose_name = "Stock Amount Unit")
    
#     cpoz_sn = models.CharField(max_length=25, blank=True, verbose_name = "CpOz SN")
#     cpoz_id = models.CharField(max_length=25, blank=True, verbose_name = "CpOz Lib ID")
#     coadd_id = models.CharField(max_length=25, blank=True, verbose_name = "CO-ADD ID")
#     chembl_id = models.CharField(max_length=25, blank=True, verbose_name = "ChEMBL ID")
#     spark_id = models.CharField(max_length=25, blank=True, verbose_name = "SPARK ID")

#     structure_id = models.ForeignKey(Chem_Structure, null=True, blank=True, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
#         db_column="structure_id", related_name="%(class)s_structureid")
#     salt_code = models.CharField(max_length=120, blank=True, verbose_name = "Salts")
#     full_mw = models.FloatField(default=0, blank=True, verbose_name = "Full MW")
#     full_mf = models.CharField(max_length=100, blank=True, verbose_name = "Full MF")
    
#     std_status = models.CharField(max_length=10, blank=True, verbose_name = "Std Status")
#     std_action = models.CharField(max_length=120, blank=True, verbose_name = "Std Action")
#     std_process = models.CharField(max_length=120, blank=True, verbose_name = "Std Process")
#     std_smiles = models.CharField(max_length=2048, blank=True, verbose_name = "Std Smiles")
 
#     pub_status = models.CharField(max_length=10, blank=True, verbose_name = "Pub Status")
#     pub_date = models.DateField(null=False, editable=False, verbose_name="Published")

# #-------------------------------------------------------------------------------------------------
# class Screen_Run(AuditModel):
#     """
#     List of Screening runs
#     """
# #-------------------------------------------------------------------------------------------------
#     HEADER_FIELDS = {
#         "run_id":"Run ID",
#         "run_type":"Run Type",
#         "assay_note":"Assa y",
#         "run_status":"Status",
#         "run_project":"Project",
#         "run_name":"Name",
#         "run_date":"Run Date",
#         "run_conditions":"Conditions",
#         "run_issues":"Issues",
#     }

#     Choice_Dictionary = {
#         'run_type':'Run_Type',
#         'run_status':'Run_Status',
#     }

#     run_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Run ID")
#     run_name = models.CharField(max_length=500, verbose_name = "Run Name")
#     run_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Type", on_delete=models.DO_NOTHING,
#         db_column="run_type", related_name="%(class)s_RunType+")
#     assay_note = models.CharField(max_length=250, blank=True, verbose_name = "Assay Note")
#     run_conditions = models.CharField(max_length=250, blank=True, verbose_name = "Run Conditions")
#     run_issues = models.CharField(max_length=250, blank=True, verbose_name = "Run Issues")
#     run_date = models.DateField(null=True, blank=True, verbose_name = "Run Date")
#     run_project = models.CharField(max_length=50, verbose_name = "Project")
#     run_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Status", on_delete=models.DO_NOTHING,
#         db_column="run_status", related_name="%(class)s_run_status+")

#     #------------------------------------------------
#     class Meta:
#         app_label = 'dscreen'
#         db_table = 'screen_run'
#         ordering=['run_type','run_id']
#         indexes = [
#             models.Index(name="run_type_idx", fields=['run_type']),
#         ]

#     # #------------------------------------------------
#     # def __str__(self) -> str:
#     #     return f"{self.run_id}"

#     #------------------------------------------------
#     def __repr__(self) -> str:
#         return f"{self.run_id} [{self.run_type}]"

#     #------------------------------------------------
#     @classmethod
#     def get(cls,RunID,verbose=0):
#         try:
#             retInstance = cls.objects.get(run_id=RunID)
#         except:
#             if verbose:
#                 print(f"[Run_ID Not Found] {RunID} ")
#             retInstance = None
#         return(retInstance)

#     #------------------------------------------------
#     @classmethod
#     def exists(cls,RunID,verbose=0):
#         return cls.objects.filter(run_id=RunID.strip()).exists()