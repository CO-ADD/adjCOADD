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

#-------------------------------------------------------------------------------------------------
# Screening Application Model
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
class Screen_Run(AuditModel):
    """
    List of Screening runs
    """
#-------------------------------------------------------------------------------------------------
    VALID_STATUS    = True
    CLASS_FIELDS = SCREENRUN_FIELDs

    Choice_Dictionary = {
        'run_type':'Run_Type',
        'screen_type':'Screen_Type',
        'run_status':'Run_Status',
    }

    run_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Run ID")
    run_name = models.CharField(max_length=500, unique=True, verbose_name = "Run Name")
    screen_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Screen Type", on_delete=models.DO_NOTHING,
        db_column="screen_type", related_name="%(class)s_ScreenType+")
    run_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Type", on_delete=models.DO_NOTHING,
        db_column="run_type", related_name="%(class)s_RunType+")
    run_conditions = models.CharField(max_length=250, blank=True, verbose_name = "Run Conditions")
    run_issues = models.CharField(max_length=250, blank=True, verbose_name = "Run Issues")
    run_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Status", on_delete=models.DO_NOTHING,
        db_column="run_status", related_name="%(class)s_RunStatus+")

    #------------------------------------------------
    class Meta:
        app_label = 'dscreen'
        db_table = 'screen_run'
        ordering=['screen_type','run_id']
        # indexes = [
        #     models.Index(name="tax_orgclass_idx", fields=['org_class']),
        #     models.Index(name="tax_taxid_idx", fields=['tax_id']),
        #     models.Index(name="tax_div_idx", fields=['division']),
        # ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.run_id} {self.run_type}"

    #------------------------------------------------
    @classmethod
    def exists(cls,RunID,verbose=0):
        try:
            retInstance = cls.objects.get(run_id=RunID)
        except:
            if verbose:
                print(f"[DataSource Not Found] {RunID} ")
            retInstance = None
        return(retInstance)
