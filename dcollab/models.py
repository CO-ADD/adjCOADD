from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.contrib.postgres.search import TrigramSimilarity
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary
from django_countries.fields import CountryField

#=================================================================================================
class Organisation(AuditModel):
    """
    List of Organisation
    """
#=================================================================================================
    LIST_VIEW_FIELDS   = {
        #"organisation_id":{'Organisation ID': {'organisation_id':URL_LINKS['organisation_id']}},
        "organisation_id":"ID",
        "organisation_code":"Code",
        "organisation_name":"Name",
        "country":"Country",
        "organisation_type":"Type",        
    }
    
    DICTIONARY_FIELDS = {
        'organisation_type':'Organisation_Type',
    }

    VIEW_GROUPS = []
    
    ID_SEQUENCE = 'Organisation'
    ID_PREFIX = 'CORG'
    ID_PAD = 5

    organisation_id = models.CharField(max_length=10, primary_key=True, verbose_name = "Organisation ID")
    organisation_code = models.CharField(max_length=20, blank=False, unique=True, verbose_name = "Organisation Code")
    organisation_name = models.CharField(max_length=250, blank=False, unique=True, verbose_name = "Organisation")
    organisation_type = models.ForeignKey(Dictionary, blank=False, verbose_name = "Organisation Type", on_delete=models.DO_NOTHING,
        db_column="organisation_type", related_name="%(class)s_organisation_type")
    country = CountryField(default='AU',verbose_name = "Country")

    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'organisation'

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.organisation_id} {self.organisation_name}"

    def __str__(self) -> str:
        return f"{self.organisation_name}"

    #------------------------------------------------
    @classmethod
    def get(cls, ID, OrganisationName=None, Code=None, verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            if ID is not None:
                retInstance = cls.objects.get(organisation_id=ID)
            elif OrganisationName is not None:
                retInstance = cls.objects.get(organisation_name=OrganisationName)
            elif Code is not None:
                retInstance = cls.objects.get(organisation_code=Code)
        except:
            if verbose:
                print(f"[Organisation Not Found] {ID} {OrganisationName} {Code} ")
            retInstance = None
        return(retInstance)

    @classmethod
    def get_bysimilarity(cls, OrganisationName=None, Similarity=0.6, verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            retInstance = cls.objects.annotate(
                            similarity=TrigramSimilarity('organisation_name', OrganisationName),
                        ).filter(similarity__gt=Similarity).order_by('-similarity').first()
        except:
            if verbose:
                print(f"[Organisation Not Found] {OrganisationName} [by Similarity] ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls, ID, OrganisationName=None, Code=None, verbose=0):
    # Returns if instance exists
        if ID is not None:
            return cls.objects.filter(organisation_id=ID).exists()
        elif OrganisationName is not None:
            return cls.objects.filter(organisation_name=OrganisationName).exists()
        elif Code is not None:
            return cls.objects.filter(organisation_code=Code).exists()
        else:
            return None
        
    # #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.organisation_id:
            self.organisation_id = self.next_id()
            if self.organisation_id: 
                super(Organisation, self).save(*args, **kwargs)
        else:
            super(Organisation, self).save(*args, **kwargs) 

#=================================================================================================
class Collab_User(AuditModel):
    """
    List of Collaborative Groups
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
        'user_id':'ID',
        'title':'Title',
        'first_name':'First Name',
        'last_name':'Last Name',
        'email':'EMail',
        'organisation_id.organisation_name':'Organisation',
        'department':'Department',
        'country.name':'Country',
    }
    DICTIONARY_FIELDS = {
    }

    ID_SEQUENCE = 'Collab_User'
    ID_PREFIX = 'CUSR'
    ID_PAD = 5

    user_id = models.CharField(max_length=15, primary_key=True, verbose_name = "User ID")
    title = models.CharField(max_length=15, blank=True, verbose_name = "Title")
    first_name = models.CharField(max_length=50, blank=True, verbose_name = "First Name")
    last_name = models.CharField(max_length=50, blank=True, verbose_name = "Last Name")
    position = models.CharField(max_length=100, blank=True, verbose_name = "Position")

    email = models.EmailField(max_length=254, blank=True, verbose_name = "EMail")
    email2 = models.EmailField(max_length=254, blank=True, verbose_name = "EMail 2nd")
    active_email = models.SmallIntegerField(default=0, blank=True, verbose_name ="Active")

    phone = models.CharField(max_length=50, blank=True, verbose_name = "Phone")
    subscribed = models.BooleanField(default=False, blank=True, verbose_name = "Newsletter")
    portal_userid = models.CharField(max_length=50, blank=True, verbose_name = "Portal UserID")
    portal_pw = models.CharField(max_length=50, blank=True, verbose_name = "Portal Password")

    organisation_id = models.ForeignKey(Organisation, null=True, blank=True, verbose_name = "Organisation ID", on_delete=models.DO_NOTHING,
        db_column="organisation_id", related_name="%(class)s_organisation_id")
        
    department = models.CharField(max_length=250, blank=True, verbose_name = "Department")
    postal_address = models.CharField(max_length=250, blank=True, verbose_name = "Postal Address")
    city = models.CharField(max_length=250, blank=True, verbose_name = "City")
    country = CountryField(default='AU',verbose_name = "Country")

    # group_id = models.ForeignKey("Collab_Group", null=True, blank=True, verbose_name = "Group Membership", on_delete=models.DO_NOTHING,
    #     db_column="group_id", related_name="%(class)s_group_id")
    pi = models.BooleanField(default=False, verbose_name = "PI")

    ora_user_id = models.CharField(max_length=15, blank=True, verbose_name = "Old User ID")
    ora_group_id = models.CharField(max_length=15, blank=True, verbose_name = "Old Group ID")
    
    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'collab_user'
        indexes = [
            models.Index(name="cuser_name_idx",fields=['first_name','last_name']),
            models.Index(name="cuser_email_idx",fields=['email']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.first_name} {self.last_name} {self.organisation_id.organisation_code}"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.first_name} {self.last_name} ({self.organisation_id.organisation_code}) [{self.user_id}]"

    #------------------------------------------------
    @classmethod
    def get(cls, ID, EMail=None, FirstName=None, LastName=None, verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            if ID is not None:
                retInstance = cls.objects.get(user_id=ID)
            elif EMail is not None:
                retInstance = cls.objects.get(email=EMail)
            elif LastName is not None:
                retInstance = cls.objects.get(first_name=FirstName, last_name=LastName)
        except:
            if verbose:
                print(f"[User Not Found] {ID} {EMail} {LastName} ")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls, ID, EMail=None, FirstName=None, LastName=None, verbose=0):
        try:
            if ID is not None:
                retValue = cls.objects.filter(user_id=ID).exists()
            elif EMail is not None:
                retValue = cls.objects.filter(email=EMail).exists()
            elif LastName is not None:
                retValue = cls.objects.filter(first_name=FirstName, last_name=LastName).exists()
        except:
            if verbose:
                print(f"[Data Not Found] {DataID} ")
            retValue = False
        return(retValue)

    # #------------------------------------------------
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
    LIST_VIEW_FIELDS = {
        'group_id':'ID',
        'group_code':'Code',
        # 'first_name':'First Name',
        # 'last_name':'Last Name',
        #'email':'E-Mail',
        'organisation_id.organisation_name':'Organisation',
        'department':'Department',
        'city':'City',
        'country.name':'Country',
        'mta_status':'MTA Status',
        'mta_document':'MTA Document'

    }
    DICTIONARY_FIELDS = {
        'mta_status':'License_Status',
    }


    VIEW_GROUPS = []

    ID_SEQUENCE = 'Collab_Group'
    ID_PREFIX = 'CGRP'
    ID_PAD = 5

    group_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Group ID")
    group_code = models.CharField(max_length=50, unique=True, verbose_name = "Group Code")

    group_members = models.ManyToManyField(Collab_User, through='Collab_Membership',through_fields=('group_id', 'user_id'))

    organisation_id = models.ForeignKey(Organisation, null=True, blank=True, verbose_name = "Organisation ID", on_delete=models.DO_NOTHING,
        db_column="organisation_id", related_name="%(class)s_organisation_id")    
    email = models.EmailField(max_length=254, blank=True, verbose_name = "EMail")
    department = models.CharField(max_length=250, blank=True, verbose_name = "Department")
    postal_address = models.CharField(max_length=250, blank=True, verbose_name = "Postal Address")
    city = models.CharField(max_length=250, blank=True, verbose_name = "City")
    country = CountryField(default='AU', verbose_name = "Country")
    pi_user_id = models.CharField(max_length=10, blank=True, verbose_name = "PI ID")
    # pi = models.ForeignKey(Collab_User, null=True, blank=True, verbose_name = "Principal Investigator", on_delete=models.DO_NOTHING,
    #     db_column="pi", related_name="%(class)s_pi")
    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_mta_status")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")

    ora_group_id = models.CharField(max_length=15, blank=True, verbose_name = "Old Group ID")
    ora_pi_id = models.CharField(max_length=15, blank=True, verbose_name = "Old PI ID")

    #------------------------------------------------
    class Meta:
        app_label = 'dcollab'
        db_table = 'collab_group'
        indexes = [
            models.Index(name="cgrp_code_idx",fields=['group_code']),
            models.Index(name="cgrp_email_idx",fields=['email']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.group_id} {self.group_code}"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.group_code} ({self.group_id})"

    #------------------------------------------------
    @classmethod
    def get(cls, ID, Code=None, PI_ID=None, Organisation_ID=None, verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            if ID is not None:
                retInstance = cls.objects.get(group_id=ID)
            elif Code is not None:
                retInstance = cls.objects.get(group_code=Code)
            elif PI_ID is not None:
                retInstance = cls.objects.get(pi_user_id=PI_ID, organisation_id=Organisation_ID)
        except:
            if verbose:
                print(f"[Group Not Found] {ID} {Code} {Organisation_ID} {PI_ID}")
            retInstance = None
        return(retInstance)

    # #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.group_id:
            self.group_id = self.next_id()
            if self.group_id: 
                super(Collab_Group, self).save(*args, **kwargs)
        else:
            super(Collab_Group, self).save(*args, **kwargs) 



#=================================================================================================
class Collab_Membership(models.Model):
    """
    List of Group Membership
    """
    MEMBERSHIP_CHOICES = [ 
            ("LI","Lead Investigator"),
            ("M","Member")
        ]

    user_id = models.ForeignKey(Collab_User, on_delete=models.CASCADE, related_name="memberships")
    group_id = models.ForeignKey(Collab_Group, on_delete=models.CASCADE, related_name="memberships")
    #date_joined = models.DateField()
    role = models.CharField(max_length=2,
            choices=MEMBERSHIP_CHOICES,
            default='M')

    class Meta:
        app_label = 'dcollab'
        db_table = 'collab_membership'
        unique_together = ('user_id', 'group_id')
        indexes = [
            models.Index(name="cmem_role_idx",fields=['role']),
        ]

    #------------------------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.group_id} <-- {self.role} -- {self.user_id}"

    def __str__(self) -> str:
        return f"{self.user_id.first_name} {self.user_id.last_name}"

    #------------------------------------------------
    @classmethod
    def get(cls, UserID, GroupID, verbose=0):
    # Returns an instance if found by ImageNAme
        try:
            retInstance = cls.objects.get(user_id=UserID, group_id=GroupID)
        except:
            if verbose:
                print(f"[Group Membership Not Found] {GroupID} <----> {UserID} ")
            retInstance = None
        return(retInstance)




#=================================================================================================
class Data_Source(AuditModel):
    """
    List of Data sources
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {

    }
    DICTIONARY_FIELDS = {
        'data_type':'Data_Type',
    }

    ID_SEQUENCE = 'Data_Source'
    ID_PREFIX = 'DSR'
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
    group_id = models.ForeignKey(Collab_Group, null=True, blank=True, verbose_name = "Group ID", on_delete=models.DO_NOTHING,
        db_column="group_id", related_name="%(class)s_group_id")

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

    # #------------------------------------------------
    # def save(self, *args, **kwargs):
    #     if not self.data_id:
    #         self.data_id = self.next_id()
    #         if self.data_id: 
    #             super(Data_Source, self).save(*args, **kwargs)
    #     else:
    #         super(Data_Source, self).save(*args, **kwargs) 


#=================================================================================================
# class Convert_UserID(AuditModel):
#     """
#     List of oraUserID -> djUserID
#     """
# #=================================================================================================

#     ora_user_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Old User ID")
#     user_id = models.CharField(max_length=15,null=False, verbose_name = "User ID")
#     project_name = models.CharField(max_length=50, blank=True, verbose_name = "Project Name")

#     class Meta:
#         app_label = 'dcollab'
#         db_table = 'convert_userid'
#         ordering=['ora_user_id']
#         indexes = [
#             models.Index(name="wusr_pid_idx", fields=['user_id']),
#             models.Index(name="wusr_pname_idx", fields=['project_name']),
# #            models.Index(name="wusr_opid_idx", fields=['old_user_id']),
#         ]


    # @classmethod
    # def new_COADD_User_ID(cls,OldUserID,verbose=0):

    #     if 'P' in OldUserID:
    #         try:
    #             _cno = int(OldProjectID[1:])
    #         except:
    #             _cno = 0
    #             return(Project.str_id(_cno))
    #         if _cno > 0 :
    #             _newID = Project.str_id(_cno)
    #             newEntry = cls()
    #             newEntry.ora_project_id = OldProjectID
    #             newEntry.project_id = _newID
    #             newEntry.save()
    #         return(newEntry)
