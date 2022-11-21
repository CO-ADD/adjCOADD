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
    
    choices_permission=(('general user', 'general user| general user can read record'), ('staff', 'staff| staff can edit record'),('admin','admin| admin is a superuser'))
#-------------------------------------------------------------------------------------------------
    username = models.CharField(primary_key=True, max_length=55,)        # J.Zuegg
    user_id = models.CharField(unique=True, max_length=50)          # uqjzuegg
    # title_name = models.CharField(max_length=15, blank=True, null=True)        # Dr
    # firstname = models.CharField(max_length=50, blank=True)        # Johannes  this field existed in AbstractUser
    # lastname = models.CharField(max_length=50, blank=True)         # Zuegg     this field existed in AbstractUser
    initials = models.CharField(max_length=5, blank=True, null=True)           # JZG
    organisation = models.CharField(max_length=250, blank=True, null=True)     # University of Queensland
    department = models.CharField(max_length=250, blank=True, null=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, blank=True, null=True)             # Blaskovich
    phone = models.CharField(max_length=25, blank=True, null=True)             # +61 7 344 62994
    # email = models.CharField(max_length=80, blank=True)             # j.zuegg@uq.edu.au    this field existed in AbstractUser
    permissions = models.CharField(max_length=250, null=True, choices=choices_permission, default=None)      # application permissions .. ReadOnly, ReadWrite, Admin ..
    session_id = models.CharField(max_length=250, blank=True, null=True)       # not sure if Django has SessionID's
    is_appuser=models.BooleanField(default=True)

    class Meta:
        db_table = 'applicationuser'
      

    def __str__(self) -> str:
        return f"{self.first_name}.{self.last_name} ({self.user_id})({self.username})"

    def save(self, *args, **kwargs):
        self.username=self.user_id
        super().save(*args, **kwargs)

#-------------------------------------------------------------------------------------------------
class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
#-------------------------------------------------------------------------------------------------
    DELETED = -9
    INVALID = -1
    UNDEFINED = 0
    VALID = 1
    CONFIRMED = 2

    astatus = models.IntegerField(verbose_name = "Status", default = 0, db_index = True, editable=False)
    acreated_at = models.DateTimeField(auto_now_add=True, null=True,verbose_name = "Created at",editable=False)
    aupdated_at = models.DateTimeField(auto_now=True, null=True,blank=True, verbose_name = "Updated at",editable=False)
    adeleted_at = models.DateTimeField(blank=True, null=True,verbose_name = "Deleted at",editable=False)
    acreated_by = models.ForeignKey(ApplicationUser, null=True, verbose_name = "Created by", related_name='%(class)s_requests_created',editable=False, on_delete=models.DO_NOTHING)
    aupdated_by = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Updated by", related_name='%(class)s_requests_updated',editable=False,on_delete=models.DO_NOTHING)
    adeleted_by = models.ForeignKey(ApplicationUser, blank=True, null=True, verbose_name = "Deleted by",related_name='%(class)s_requests_deleted',editable=False,on_delete=models.DO_NOTHING)

    class Meta:
        abstract = True
    
    def delete(self,**kwargs):
        print('this is delete function from Audit model class')
        self.astatus = -9
        self.adeleted_at = timezone.now()
        self.adeleted_by = kwargs.get("user")
        print(f"deleted by {self.adeleted_by}")
        self.save(**kwargs)

    def save(self, *args, **kwargs):
        print('this is save function from Audit model class')
        user=kwargs.get("user")
        print(user)
        if not self.acreated_by: 	#Createing
            self.acreated_by = user       
        else:					#Updateing
            self.aupdated_by = user
        if kwargs.get("user"):
            kwargs.pop("user")
        super().save(*args, **kwargs)

#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
class Dictionary(AuditModel):
#-------------------------------------------------------------------------------------------------
    
    Dictionary_Class= models.CharField(max_length=30, verbose_name = "Dictionary_Class")
    Dict_Value =models.CharField(primary_key=True, unique=True, max_length=50, verbose_name = "Value"  )
    Dict_Desc = models.CharField(max_length=120, blank=True, null=True, verbose_name = "Description")
   
    def __str__(self) -> str:
        return f"{self.Dict_Value}.{self.Dict_Desc}"
    
    class Meta:
        db_table = 'dictionary'


#-------------------------------------------------------------------------------------------------
class ApplicationLog(models.Model):
#-------------------------------------------------------------------------------------------------
    log_code = models.CharField(max_length=15, blank=True, editable=False)
    log_proc = models.CharField(max_length=50, blank=True, editable=False)
    log_type = models.CharField(max_length=15, blank=True, editable=False)
    log_time = models.DateTimeField(auto_now=True, blank=True, editable=False)
    log_user = models.ForeignKey(ApplicationUser, editable=False, on_delete=models.DO_NOTHING)
    log_object = models.CharField(max_length=15, blank=True, editable=False)
    log_desc = models.CharField(max_length=1024, blank=True, editable=False)
    log_status = models.CharField(max_length=15, blank=True, editable=False)

    class Meta:
        db_table = 'applicationLog'



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