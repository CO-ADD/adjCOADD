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

#=================================================================================================
class Organisation(AuditModel):
    """
    List of Organisation
    """
#=================================================================================================
    HEADER_FIELDS   = {}
    Choice_Dictionary = {
        'org_type':'Organisation_Type',
    }

    organisation = models.CharField(max_length=250, blank=True, verbose_name = "Organisation")
    org_code = models.CharField(max_length=10, blank=True, verbose_name = "Organisation Code")
    org_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Organisation Type", on_delete=models.DO_NOTHING,
        db_column="org_type", related_name="%(class)s_OrgType+")
    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'organisation'

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.organisation}"

#=================================================================================================
class Collab_User(AuditModel):
    """
    List of Collaborative Groups
    """
#=================================================================================================
    HEADER_FIELDS = {}
    Choice_Dictionary = {}

    user_id = models.CharField(max_length=10, primary_key=True, verbose_name = "User ID")
    title = models.CharField(max_length=15, blank=True, verbose_name = "Title")
    first_name = models.CharField(max_length=50, blank=True, verbose_name = "First Code")
    last_name = models.CharField(max_length=50, blank=True, verbose_name = "Last Code")
    position = models.CharField(max_length=100, blank=True, verbose_name = "Position")
    email = models.CharField(max_length=50, blank=True, verbose_name = "EMail")
    phone = models.CharField(max_length=50, blank=True, verbose_name = "Phone")
    subscribed = models.BooleanField(default=False, blank=True, verbose_name = "Newsletter")
    portal_userid = models.CharField(max_length=50, blank=True, verbose_name = "UserID")
    portal_pw = models.CharField(max_length=50, blank=True, verbose_name = "Password")
    organisation = models.ForeignKey(Organisation, null=True, blank=True, verbose_name = "Organisation", on_delete=models.DO_NOTHING,
        db_column="organisation", related_name="%(class)s_Organisation+")    
    department = models.CharField(max_length=250, blank=True, verbose_name = "Department")
    postal_address = models.CharField(max_length=250, blank=True, verbose_name = "Postal Address")
    city = models.CharField(max_length=250, blank=True, verbose_name = "City")
    country = models.CharField(max_length=250, blank=True, verbose_name = "Country")
    group = models.ForeignKey("Collab_Group", null=True, blank=True, verbose_name = "Group Membership", on_delete=models.DO_NOTHING,
        db_column="group", related_name="%(class)s_Group+")
    
    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'collab_user'

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.first_name} {self.last_name} {self.organisation.org_code}"

#=================================================================================================
class Collab_Group(AuditModel):
    """
    List of Collaborative Groups
    """
#=================================================================================================
    HEADER_FIELDS = {}
    Choice_Dictionary = {
        'mta_status':'License_Status',
    }

    group_id = models.CharField(max_length=10, primary_key=True, verbose_name = "Group ID")
    group_code = models.CharField(max_length=10, unique=True, verbose_name = "Group Code")
    organisation = models.ForeignKey(Organisation, null=True, blank=True, verbose_name = "Organisation", on_delete=models.DO_NOTHING,
        db_column="organisation", related_name="%(class)s_Organisation+")
    department = models.CharField(max_length=250, blank=True, verbose_name = "Department")
    postal_address = models.CharField(max_length=250, blank=True, verbose_name = "Postal Address")
    city = models.CharField(max_length=250, blank=True, verbose_name = "City")
    country = models.CharField(max_length=250, blank=True, verbose_name = "Country")
    pi = models.ForeignKey(Collab_User, null=True, blank=True, verbose_name = "Principal Investigator", on_delete=models.DO_NOTHING,
        db_column="pi", related_name="%(class)s_PI+")
    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_MTA+")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'collab_group'

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.group_code}"

#=================================================================================================
class Data_Source(AuditModel):
    """
    List of Data sources
    """
#=================================================================================================
    HEADER_FIELDS = {}
    Choice_Dictionary = {
        'source_type':'DataSource_Type',
    }

    source_id = models.CharField(max_length=25,primary_key=True, verbose_name = "Source ID")
    source_name = models.CharField(max_length=50, blank=True, verbose_name = "Source Name")
    source_code = models.CharField(max_length=10, blank=True, verbose_name = "Source Code")
    description = models.CharField(max_length=1000, blank=True, verbose_name = "Description")
    source_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Source Type", on_delete=models.DO_NOTHING,
        db_column="source_type", related_name="%(class)s_SourceType+")
    journal = models.CharField(max_length=50, blank=True, verbose_name = "Journal")
    year = models.IntegerField(blank=True, verbose_name = "Year")
    volume = models.CharField(max_length=50, blank=True, verbose_name = "Volume")
    issue = models.CharField(max_length=50, blank=True, verbose_name = "Issue")
    page = models.CharField(max_length=50, blank=True, verbose_name = "Page")
    title = models.CharField(max_length=500, blank=True, verbose_name = "Title")
    pubmed_id = models.IntegerField(blank=True, verbose_name = "PubMed ID")
    doi = models.CharField(max_length=100, blank=True, verbose_name = "DOI")
    url = models.CharField(max_length=100, blank=True, verbose_name = "URL")
    authors = models.CharField(max_length=1000, blank=True, verbose_name = "Authors")
    patent_id = models.CharField(max_length=500, blank=True, verbose_name = "Patent IDs")
    collab_group = models.ForeignKey(Collab_Group, null=True, blank=True, verbose_name = "Group", on_delete=models.DO_NOTHING,
        db_column="collab_group", related_name="%(class)s_CollabGroup+")

    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'data_source'
        ordering=['source_name','source_type']
        # indexes = [
        #     models.Index(name="tax_orgclass_idx", fields=['org_class']),
        #     models.Index(name="tax_taxid_idx", fields=['tax_id']),
        #     models.Index(name="tax_div_idx", fields=['division']),
        # ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.source_id} {self.source_type} {self.source_name}"

    #------------------------------------------------
    @classmethod
    def exists(cls,SourceID,verbose=0):
        try:
            retInstance = cls.objects.get(source_id=SourceID)
        except:
            if verbose:
                print(f"[DataSource Not Found] {SourceID} ")
            retInstance = None
        return(retInstance)
