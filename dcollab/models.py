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

    ID_SEQUENCE = 'Organisation'
    ID_PREFIC = 'ORG'
    ID_PAD = 5

    org_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Organisation ID")
    org_name = models.CharField(max_length=250, blank=True, verbose_name = "Organisation")
    org_code = models.CharField(max_length=10, blank=True, verbose_name = "Organisation Code")
    org_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Organisation Type", on_delete=models.DO_NOTHING,
        db_column="org_type", related_name="%(class)s_OrgType+")
    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'organisation'

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.org_id} {self.org_name}"

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.org_id:
            self.org_id = self.next_id()
            if self.org_id: 
                super(Organisation, self).save(*args, **kwargs)
        else:
            super(Organisation, self).save(*args, **kwargs) 

#=================================================================================================
class Collab_User(AuditModel):
    """
    List of Collaborative Groups
    """
#=================================================================================================
    HEADER_FIELDS = {}
    Choice_Dictionary = {}

    ID_SEQUENCE = 'Collab_User'
    ID_PREFIC = 'CUSR'
    ID_PAD = 5

    user_id = models.CharField(max_length=15, primary_key=True, verbose_name = "User ID")
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
    def __repr__(self) -> str:
        return f"{self.first_name} {self.last_name} {self.organisation.org_code}"

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.user_id:
            self.user_id = self.next_id()
            if self.user_id: 
                super(Collab_User, self).save(*args, **kwargs)
        else:
            super(Collab_User, self).save(*args, **kwargs) 

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

    ID_SEQUENCE = 'Collab_Group'
    ID_PREFIC = 'CGRP'
    ID_PAD = 5

    group_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Group ID")
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
    def __repr__(self) -> str:
        return f"{self.group_code}"

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.group_id:
            self.group_id = self.next_id()
            if self.group_id: 
                super(Collab_Group, self).save(*args, **kwargs)
        else:
            super(Collab_Group, self).save(*args, **kwargs) 

#=================================================================================================
class Data_Source(AuditModel):
    """
    List of Data sources
    """
#=================================================================================================
    HEADER_FIELDS = {

    }
    Choice_Dictionary = {
        'data_type':'Data_Type',
    }

    ID_SEQUENCE = 'Data_Source'
    ID_PREFIC = 'DSR'
    ID_PAD = 5

    data_id = models.CharField(max_length=25,primary_key=True, verbose_name = "Data ID")
    data_name = models.CharField(max_length=50, blank=True, verbose_name = "Data Name")
    data_code = models.CharField(max_length=10, blank=True, verbose_name = "Data Code")
    description = models.CharField(max_length=1000, blank=True, verbose_name = "Description")
    data_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Data Type", on_delete=models.DO_NOTHING,
        db_column="data_type", related_name="%(class)s_data_type")
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
        ordering=['data_name','data_type']
      
    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.data_id} "

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.data_id} {self.data_code} {self.data_type}"

    #------------------------------------------------
    @classmethod
    def exists(cls,DataID,verbose=0):
        try:
            retInstance = cls.objects.get(data_id=DataID)
        except:
            if verbose:
                print(f"[Data Not Found] {DataID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.data_id:
            self.data_id = self.next_id()
            if self.data_id: 
                super(Data_Source, self).save(*args, **kwargs)
        else:
            super(Data_Source, self).save(*args, **kwargs) 
