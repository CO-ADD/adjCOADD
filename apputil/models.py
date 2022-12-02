# Create your models here.
from django.db import models
from django.contrib.auth.models import AbstractUser, Group, Permission
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse
from django import forms
from django.utils import timezone
import logging
# Create your models here.


#-------------------------------------------------------------------------------------------------
class ApplicationUser(AbstractUser):    
#-------------------------------------------------------------------------------------------------
    username = models.CharField(unique=True, max_length=55, verbose_name='uquser')       # uqjzuegg 
    name = models.CharField(primary_key=True,  max_length=50, verbose_name='user')          # J.Zuegg
    initials = models.CharField(max_length=5, null=True, blank=True)           # JZG
    organisation = models.CharField(max_length=250, null=True, blank=True)     # University of Queensland
    department = models.CharField(max_length=250, null=True, blank=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, null=True, blank=True)             # Blaskovich
    phone = models.CharField(max_length=25, null=True, blank=True)             # +61 7 344 62994
    permission = models.CharField(max_length=10, default = 'No', null=False)      # application permissions .. Read, Write, Delete, Admin ..
    is_appuser=models.BooleanField(default=True)

    class Meta:
        db_table = 'app_user'
    
    # --------------------------------------------------------------------------
    def has_permission(self,strPermission) -> bool:
    #
    # Returns True/False if User has strPermission 
    #   checks if strPermission/self.permission in [Read, Write, Delete, Admin]
    # --------------------------------------------------------------------------
        _Permissions = {
            'Read':1,
            'Write':2,
            'Delete':3,
            'Admin':10,
        }
        if strPermission in _Permissions:
            if self.permission in _Permissions: 
                return(_Permissions[self.permission]>=_Permissions[strPermission])
            else:
                return(False)
        else:
            return(False)

    # --------------------------------------------------------------------------
    def __str__(self) -> str:
        return f"{self.name}" 


#-------------------------------------------------------------------------------------------------
class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
#-------------------------------------------------------------------------------------------------
    DELETED   = -9
    INVALID   = -1
    UNDEFINED = 0
    VALID     = 1
    CONFIRMED = 2
    OWNER     = "orgdb"

    astatus = models.IntegerField(verbose_name = "Status", default = 0, db_index = True, editable=False)
    acreated_at = models.DateTimeField(null=False, editable=False, verbose_name="Created at")
    aupdated_at = models.DateTimeField(null=True,  editable=False, verbose_name="Updated at")
    adeleted_at = models.DateTimeField(null=True,  editable=False, verbose_name="Deleted at",)
    acreated = models.ForeignKey(ApplicationUser, null=False, verbose_name = "Created by", 
        related_name="%(class)s_acreated_by", editable=False, on_delete=models.DO_NOTHING)
    aupdated = models.ForeignKey(ApplicationUser, null=True,  verbose_name = "Updated by", 
        related_name="%(class)s_aupdated_by", editable=False, on_delete=models.DO_NOTHING)
    adeleted = models.ForeignKey(ApplicationUser, null=True,  verbose_name = "Deleted by", 
        related_name="%(class)s_adeleted_by", editable=False, on_delete=models.DO_NOTHING)

       #------------------------------------------------
    class Meta:
        abstract = True
    
    #------------------------------------------------
    def delete(self,**kwargs):
        appuser=kwargs.get("user")
        kwargs.pop("user")
        if appuser is None:
            appuser = ApplicationUser.objects.get(name=self.OWNER)

        self.astatus = self.DELETED
        self.adeleted_id = appuser
        self.adeleted_at = timezone.now()
        super(AuditModel,self).save(**kwargs)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        appuser=kwargs.get("user")
        kwargs.pop("user")
        if appuser is None:
            appuser = ApplicationUser.objects.get(name=self.OWNER)

        if not self.acreated_id:
            self.acreated_id = appuser
            self.acreated_at = timezone.now()       
        else:	
            self.aupdated_id = appuser
            self.aupdated_at = timezone.now()       
        super(AuditModel,self).save(*args, **kwargs)

#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
class Dictionary(AuditModel):
#-------------------------------------------------------------------------------------------------
    
    dict_value =models.CharField(primary_key=True, unique=True, max_length=50, verbose_name = "Value"  )
    dict_class= models.CharField(max_length=30, verbose_name = "Class")
    dict_desc = models.CharField(max_length=150, blank=True, null=True, verbose_name = "Description")
   
    #------------------------------------------------
    class Meta:
        db_table = 'dictionary'
        ordering=['dict_value']
        indexes = [
            models.Index(name="dict_class_idx",fields=['dict_class']),
        ]
    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.dict_value}.{self.dict_desc}"


#-------------------------------------------------------------------------------------------------
class ApplicationLog(models.Model):
#-------------------------------------------------------------------------------------------------
    log_code = models.CharField(max_length=15, null=True, blank=True, editable=False)
    log_proc = models.CharField(max_length=50, null=True, blank=True, editable=False)
    log_type = models.CharField(max_length=15, null=True, blank=True, editable=False)
    log_time = models.DateTimeField(auto_now=True, editable=False)
    log_user = models.ForeignKey(ApplicationUser, db_column = "log_user", editable=False, on_delete=models.DO_NOTHING)
    log_object = models.CharField(max_length=15, null=True, blank=True, editable=False)
    log_desc = models.CharField(max_length=1024, null=True, blank=True, editable=False)
    log_status = models.CharField(max_length=15, null=True, blank=True, editable=False)

    class Meta:
        db_table = 'app_log'



#=========================Not used...==============
class ChoiceArrayField(ArrayField):
    """
    A field that allows us to store an array of choices.
    
    Uses Django 1.9's postgres ArrayField
    and a MultipleChoiceField for its formfield.
    
    Usage:
        
        choices = ChoiceArrayField(models.CharField(max_length=...,
                                                    choices=(...,)),
                                   default=[...])
    """
#-------------------------------------------------------------------------------------------------
    def formfield(self, **kwargs):
        defaults = {
            'form_class': forms.MultipleChoiceField,
            'choices': self.base_field.choices,
        }
        defaults.update(kwargs)
  
        return super(ArrayField, self).formfield(**defaults)

class Integer(models.Model):
    LOW = 5
    NORMAL = 10
    HIGH = 25
    STATUS_CHOICES = (
        (LOW, 'Low'),
        (NORMAL, 'Normal'),
        (HIGH, 'High'), 
    )
    num= models.IntegerField(default = 5, choices=STATUS_CHOICES)