from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *


from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary

#=================================================================================================
class Project(AuditModel):
    """
    List of Projects, Chem Library 
    """
#=================================================================================================
    Choice_Dictionary = {
        'project_class':'Project_Class',
    }
    ID_SEQUENCE = 'Project'
    ID_PREFIX = 'P'
    ID_PAD = 5

    project_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Project ID")
    project_name = models.CharField(max_length=50, unique=True, verbose_name = "Project Name")
    project_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="project_class", related_name="%(class)s_projectclass")
    project_owner =  models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="project_class", related_name="%(class)s_projectclass")
    
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    class Meta:
        app_label = 'dscreen'
        db_table = 'project'
        ordering=['project_id']
        indexes = [
            models.Index(name="prj_pname_idx", fields=['project_name']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.project_id}  {self.source}"

    #------------------------------------------------
    @classmethod
    def get(cls,ProjectID,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            retInstance = cls.objects.get(project_id=ProjectID)
        except:
            retInstance = None
            if verbose:
                print(f"[Project Not Found] {ProjectID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,ProjectID,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        retValue = cls.objects.filter(project_id=ProjectID).exists()
        return(retValue)


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.project_id:
            self.project_id = self.next_id()
            if self.project_id: 
                super(Project, self).save(*args, **kwargs)
        else:
            super(Project, self).save(*args, **kwargs) 


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

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.run_id}"

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

