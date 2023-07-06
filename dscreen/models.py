from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *


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
    HEADER_FIELDS = {
        "run_id":"RunID",
        "run_type":"RunType",
        "run_status":"Status",
        "run_owner":"Owner",
        "run_name":"Name",
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
    run_conditions = models.CharField(max_length=250, blank=True, verbose_name = "Run Conditions")
    run_issues = models.CharField(max_length=250, blank=True, verbose_name = "Run Issues")
    run_owner = models.CharField(max_length=500, verbose_name = "Run Owner")
    run_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Run Status", on_delete=models.DO_NOTHING,
        db_column="run_status", related_name="%(class)s_RunStatus+")

    #------------------------------------------------
    class Meta:
        app_label = 'dscreen'
        db_table = 'screen_run'
        ordering=['screen_type','run_id']
        indexes = [
            models.Index(name="run_type_idx", fields=['run_type']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.run_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.run_id} [{self.run_type}]"

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

    #------------------------------------------------
    @classmethod
    def exists(cls,RunID,verbose=0):
        return cls.objects.filter(run_id=RunID.strip()).exists()

