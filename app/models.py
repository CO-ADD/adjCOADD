
# Create your models here.
from django.db import models
from django.contrib.auth.models import AbstractUser, Group, Permission
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse
from django import forms
# Create your models here.


#-------------------------------------------------------------------------------------------------
class ApplicationUser(AbstractUser):
#-------------------------------------------------------------------------------------------------
    user_id = models.CharField(unique=True, max_length=50)          # uqjzuegg
    title_name = models.CharField(max_length=15, blank=True)        # Dr
    # firstname = models.CharField(max_length=50, blank=True)        # Johannes  this field existed in AbstractUser
    # lastname = models.CharField(max_length=50, blank=True)         # Zuegg     this field existed in AbstractUser
    short_name = models.CharField(max_length=55, blank=True)        # J.Zuegg
    initials = models.CharField(max_length=5, blank=True)           # JZG
    organisation = models.CharField(max_length=250, blank=True)     # University of Queensland
    department = models.CharField(max_length=250, blank=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, blank=True,)             # Blaskovich
    phone = models.CharField(max_length=25, blank=True)             # +61 7 344 62994
    # email = models.CharField(max_length=80, blank=True)             # j.zuegg@uq.edu.au    this field existed in AbstractUser
    permissions = models.CharField(max_length=250, blank=True)      # application permissions .. ReadOnly, ReadWrite, Admin ..
    session_id = models.CharField(max_length=250, blank=True)       # not sure if Django has SessionID's
    is_appuser=models.BooleanField(default=True)

    class Meta:
        db_table = 'applicationuser'
        # permissions = (('update', 'change data'),)

    def __str__(self) -> str:
        return f"{self.first_name}.{self.last_name} ({self.user_id})"


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
        self.astatus = -9
        # self.adeleted_at = timezone.now()
        self.adeleted_by = kwargs.get("user")
        self.save(updated_fields = ['adeleted_at','adeleted_by','astatus'])
    
    def save(self, *args, **kwargs):
        user = kwargs.get("user")
        if self.pk: #Object already exists
            self.aupdated_by = user
        else:
            self.acreated_by = user
        super(AuditModel,self).save(*args, **kwargs)


#-------------------------------------------------------------------------------------------------
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
        # Skip our parent's formfield implementation completely as we don't
        # care for it.
        # pylint:disable=bad-super-call
        return super(ArrayField, self).formfield(**defaults)


#-------------------------------------------------------------------------------------------------
class Dictionaries(AuditModel):
#-------------------------------------------------------------------------------------------------


    Dictionary_ID = models.CharField(max_length=30, db_index=True, verbose_name = "Dictionary", )
    Dictionary_Class= models.CharField(max_length=30, verbose_name = "Dictionary_Class", )
    Dict_Value =models.CharField(max_length=50, verbose_name = "Value",  ) # ArrayField(models.CharField(max_length=150, null=True, blank=True,),  default=list, blank=True)  #
    Dict_Desc = models.CharField(max_length=120, blank=True, null=True, verbose_name = "Description")
    Dict_Value_Type = models.CharField(max_length=20, verbose_name = "Type")
    Dict_View_Order = models.IntegerField(verbose_name = "View Order", null=True)

    def __str__(self) -> str:
        return f"{self.Dict_Value}.{self.Dict_Desc}"


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


class Mytest2(models.Model):
    
    tag=models.CharField(max_length=155, default='default') 
    date=models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.tag
