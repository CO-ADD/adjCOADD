
# Create your models here.
from django.db import models
from django.contrib.auth.models import AbstractUser, Group, Permission
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse

# Create your models here.



class User(AbstractUser):
    role=models.CharField(max_length=250, null=True)

class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
    astatus = models.IntegerField(verbose_name = "Status", default = 0,  editable=False, db_index=True)   #index = True,
    acreated_at = models.DateTimeField(auto_now_add=True, null=True, verbose_name = "Created at",editable=False)
    aupdated_at = models.DateTimeField(auto_now=True, null=True, verbose_name = "Updated at",editable=False)
    adeleted_at = models.DateTimeField(null=True, verbose_name = "Deleted at",editable=False)
    acreated_by = models.ForeignKey(User,blank=True, null=True, verbose_name = "Created by", related_name='%(class)s_requests_created', editable=False, on_delete=models.CASCADE) #%(class)s_
    aupdated_by = models.ForeignKey(User, blank=True,  null=True, verbose_name = "Updated by", related_name='%(class)s_requests_updated',editable=False, on_delete=models.CASCADE)
    adeleted_by = models.ForeignKey(User, blank=True, null=True, verbose_name = "Deleted by",related_name='%(class)s_requests_deleted',editable=False, on_delete=models.CASCADE)

    class Meta:
        abstract = True
        indexes = [
            models.Index(fields=['astatus']),
            
        ]
    
    def delete(self,**kwargs):
        self.astatus = -9
        self.adeleted_at = timezone.now()
        self.adeleted_by = kwargs.get("user")
        self.save(updated_fields = ['adeleted_at','adeleted_by','astatus'])
    
    def save(self, *args, **kwargs):
        user = kwargs.get("user")
        if self.pk: #Object already exists
            self.aupdated_by = user
        else:
            self.acreated_by = user
        super(AuditModel,self).save(*args, **kwargs)

class ApplicationDictionary(AuditModel):
    dict_table = models.CharField(max_length=30, verbose_name = "Dict Table")
    dict_field = models.CharField(max_length=15,  verbose_name = "Dict Field")
    dict_value = models.CharField(max_length=50, verbose_name = "Dict Value")
    dict_value_type = models.CharField(max_length=20, verbose_name = "Dict Value Type")
    dict_order = models.IntegerField(verbose_name = "Dict Value Order", null=True, blank=True)
    dict_desc = models.CharField(max_length=120, blank=True, verbose_name = "Dict Value Description")

    class Meta:
        db_table = 'application_dictionary'
        unique_together = (('dict_value_type', 'dict_value'),)

    def __str__(self) -> str:
        return f"{self.dict_value} ({self.dict_table}.{self.dict_field})"


class ApplicationLog(models.Model):
    log_code = models.CharField(max_length=15, blank=True, editable=False)
    log_proc = models.CharField(max_length=50, blank=True, editable=False)
    log_type = models.CharField(max_length=15, blank=True, editable=False)
    log_time = models.DateTimeField(auto_now=True, blank=True, editable=False)
    log_user = models.ForeignKey(User, editable=False, on_delete=models.CASCADE)
    log_object = models.CharField(max_length=15, blank=True, editable=False)
    log_desc = models.CharField(max_length=1024, blank=True, editable=False)
    log_status = models.CharField(max_length=15, blank=True, editable=False)

    class Meta:
        db_table = 'application_log'


class ApplicationUser(AuditModel):
    user_id = models.CharField(unique=True, max_length=50)          # uqjzuegg
    title_name = models.CharField(max_length=15, blank=True)        # Dr
    first_name = models.CharField(max_length=50, blank=True)        # Johannes
    last_name = models.CharField(max_length=50, blank=True)         # Zuegg
    short_name = models.CharField(max_length=55, blank=True)        # J.Zuegg
    initials = models.CharField(max_length=5, blank=True)           # JZG
    organisation = models.CharField(max_length=250, blank=True)     # University of Queensland
    department = models.CharField(max_length=250, blank=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, blank=True)             # Blaskovich
    phone = models.CharField(max_length=25, blank=True)             # +61 7 344 62994
    email = models.CharField(max_length=80, blank=True)             # j.zuegg@uq.edu.au
    permissions = models.CharField(max_length=250, blank=True)      # application permissions .. ReadOnly, ReadWrite, Admin ..
    session_id = models.CharField(max_length=250, blank=True)       # not sure if Django has SessionID's

    class Meta:
        db_table = 'application_user'
        permissions = (('update', 'change data'),)

    def __str__(self) -> str:
        return f"{self.first_name}.{self.last_name} ({self.user_id})"
  
