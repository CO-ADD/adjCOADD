import pandas as pd
import numpy as np
from asgiref.sync import sync_to_async

from django.db import models
from django.core.exceptions import ValidationError
from django.contrib.auth.models import AbstractUser, Group, Permission
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse
from django import forms
from django.utils import timezone
import logging

from adjcoadd.constants import *
# Create your models here.


#-------------------------------------------------------------------------------------------------
class ApplicationUser(AbstractUser):    
#-------------------------------------------------------------------------------------------------
    username = models.CharField(unique=True, max_length=55, verbose_name='uq user')       # uqjzuegg 
    name = models.CharField(primary_key=True,  max_length=50, verbose_name='user name')          # J.Zuegg
    initials = models.CharField(max_length=5, null=True, blank=True)           # JZG
    organisation = models.CharField(max_length=250, null=True, blank=True)     # University of Queensland
    department = models.CharField(max_length=250, null=True, blank=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, null=True, blank=True)             # Blaskovich
    phone = models.CharField(max_length=25, null=True, blank=True)             # +61 7 344 62994
    permission = models.CharField(max_length=10, default = 'No', null=False)      # application permissions .. Read, Write, Delete, Admin ..
    is_appuser=models.BooleanField(default=True)

    class Meta:
        db_table = 'app_user'
        ordering=['username']
    
    #------------------------------------------------
    @classmethod
    #
    # Returns an User instance if found by name
    #
    def exists(self,UserName):
        try:
            retInstance = self.objects.get(name=UserName)
        except:
            retInstance = None
        return(retInstance)

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
    
    # get field names in postgres in the order provided by constants.py
    @classmethod
    def get_databasefields(self, fields=None):
        if fields:
            databasefields=fields.keys()
        else:
            databasefields=None
        return databasefields
    # get field verbose or customized name in the order provided by constants.py
    @classmethod
    def get_fields(self, fields=None):
        if fields:
            select_fields=[fields[f.name] for f in self._meta.fields if f.name in fields.keys()]
        else:
            select_fields=None
        return select_fields


    # get field name in model Class in the order provided by constants.py
    @classmethod
    def get_modelfields(self, fields=None):
        if fields:
            model_fields=[f.name for f in self._meta.fields if f.name in fields.keys()]
        else:
            model_fields=None
        return model_fields
    #------------------------------------------------
    
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
    UNDEFINED =  0
    VALID     =  1
    CONFIRMED =  2
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
    def validate(self,**kwargs):
    #
    # Validates the instance using full_clean
    # 
        retValid = None
        try:
            self.full_clean(**kwargs)
        except ValidationError as e:
            retValid = e.message_dict
        return(retValid)


    #-------------------------------------------------------------------
    def clean_Fields(self,default_Char="",default_Integer=0):
    #
    # Sets 'None' fields in the instance according to Django guidelines 
    #   sets CharField    to "" (empty) or 'default' 
    #   sets IntegerField to 0 or 'default'
    #
        clFields = {}
        for field in self._meta.get_fields(include_parents=False):
            fType = field.get_internal_type()
            if fType == "IntegerField":
                if hasattr(self,field.name):
                    if getattr(self,field.name) is None:
                        defValue = default_Integer
                        fDict = field.deconstruct()[3]
                        if 'default' in fDict:
                            defValue = fDict['default']
                        setattr(self,field.name,defValue)
                        clFields[field.name]=defValue
            elif fType == "CharField":
                if hasattr(self,field.name):
                    if getattr(self,field.name) is None:
                        defValue = default_Char
                        fDict = field.deconstruct()[3]
                        if 'default' in fDict:
                            defValue = fDict['default']
                        setattr(self,field.name,defValue)
                        clFields[field.name]=defValue
        return(clFields)

    #------------------------------------------------
    def delete(self,**kwargs):
        appuser=kwargs.get("user")
        kwargs.pop("user",None)
        if appuser is None:
            appuser = ApplicationUser.objects.get(name=self.OWNER)

        self.astatus = self.DELETED
        print(f'delete item {self.astatus}')
        self.adeleted_id = appuser
        self.adeleted_at = timezone.now()
        super(AuditModel,self).save(**kwargs)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        #
        # Checks for application user
        #   could also use middleware
        #
        appuser=kwargs.get("user")
        kwargs.pop("user",None)
        if appuser is None:
            appuser = ApplicationUser.objects.get(name=self.OWNER)

        if not self.acreated_id:
            self.acreated_id = appuser
            self.acreated_at = timezone.now()       
        else:	
            self.aupdated_id = appuser
            self.aupdated_at = timezone.now()

        #
        # Checks if a clean=True is requested
        #   Default, via forms is False, but can be set via scripts/API
        #  
        modelClean=kwargs.get("clean")  
        kwargs.pop("clean",None)
        if modelClean:
            self.full_clean()
                 
        super(AuditModel,self).save(*args, **kwargs)

    #------------------------------------------------
    #Method Get Fields, Values List
    # get field names in postgres in the order provided by constants.py
    @classmethod
    def get_databasefields(self, fields=None):
        if fields:
            databasefields=fields.keys()
        else:
            databasefields=None
        return databasefields
    # get field verbose or customized name in the order provided by constants.py
    @classmethod
    def get_fields(self, fields=None):
        if fields:
            select_fields=[fields[f.name] for f in self._meta.fields if f.name in fields.keys()]
        else:
            select_fields=None
        return select_fields
    #------------------------------------------------
    # get field name in model Class in the order provided by constants.py
    @classmethod
    def get_modelfields(self, fields=None):
        if fields:
            model_fields=[f.name for f in self._meta.fields if f.name in fields.keys()]
        else:
            model_fields=None
        return model_fields
 
    #------------------------------------------------
    # objects values according to fields return from the above class methods
    def get_values(self, fields=None):
        value_list=[]
        for field in self._meta.fields:
            if field.name in fields.keys():
                obj=getattr(self, field.name)
                if obj:
                    if isinstance(obj, list):
                        array_to_string=','.join(str(e) for e in obj)
                        value_list.append(array_to_string)
                    else:   
                        value_list.append(field.value_to_string(self))
                else:
                    value_list.append(" ")
        return value_list
#-------------------------------------------------------------------------------------------------
    # data-visulization
    @classmethod
    # @sync_to_async
    def get_pivottable(self, querydata, columns_str, index_str,aggfunc, values):
        np_aggfunc={"Sum": np.sum, "Mean":np.mean, "Std":np.std}
        data=list(querydata.values())
        df=pd.DataFrame(data)
        columns=columns_str.split(",") 
        index=index_str.split(",")
        table=pd.pivot_table(df, values=values, index=index,
                        columns=columns, aggfunc=np_aggfunc[aggfunc])
        return table



#-------------------------------------------------------------------------------------------------
class Dictionary(AuditModel):
#-------------------------------------------------------------------------------------------------
    
    dict_value =models.CharField(primary_key=True, unique=True, max_length=50, verbose_name = "Value"  )
    dict_class= models.CharField(max_length=30, verbose_name = "Class")
    dict_desc = models.CharField(max_length=120, blank=True, verbose_name = "Description")
   
    #------------------------------------------------
    class Meta:
        app_label = 'apputil'
        db_table = 'app_dictionary'
        ordering=['dict_class','dict_value']
        indexes = [
            models.Index(name="dict_class_idx",fields=['dict_class']),
        ]
    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.dict_value} <small class='not-visible'>({self.dict_desc})</small>"


    #------------------------------------------------
    @classmethod
    def exists(self,DictClass,DictValue=None,DictDesc=None,verbose=1):
    #
    # Returns a Dictionary instance if found 
    #    by dict_value
    #    by dict_desc (set dict_value = None)
    #
        if DictValue:
            try:
                retDict = Dictionary.objects.get(dict_value=DictValue, dict_class=DictClass)
            except:
                if verbose:
                    print(f"[Dict Value Not Found] {DictValue} {DictClass}")
                retDict = None
        elif DictDesc:
            try:
                retDict = Dictionary.objects.get(dict_desc=DictDesc, dict_class=DictClass)
            except:
                if verbose:
                    print(f"[Dict Desc Not Found] {DictDesc} {DictClass}")
                retDict = None
        else:
            retDict = None
        return(retDict)

    #------------------------------------------------
    @classmethod
    #
    # Returns Dictionary entries for a DictClass as Choices
    #
    def get_Dictionary_asChoice(self, DictClass, sep = " | ", emptyChoice= ('--', 'empty')):
        dictList=Dictionary.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc')
        if dictList:
            choices_test=tuple([tuple(d.values()) for d in dictList])
            choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_test)    
        else:
            choices=emptyChoice
        return choices

    #------------------------------------------------
    @classmethod
    #
    # Returns Dictionary entries for a DictClass as Choices
    #
    def get_DictionaryStrList_asArray(self,DictClass,DictValueStr=None,DictDescStr=None,sep=";",notFound="#"):
    #-----------------------------------------------------------------------------------
        retDictList = []
        if DictValueStr:
            dLst = DictValueStr.split(sep)
            for dVal in dLst:
                xDict = self.exists(DictClass,dVal.strip(),None)
                if xDict:
                    retDictList.append(xDict.dict_value)
                else:
                    retDictList.append(f"{dVal.strip()}{notFound}")
        elif DictDescStr:
            dLst = DictDescStr.split(sep)
            for dDesc in dLst:
                xDict = self.exists(DictClass,None,dDesc.strip())
                if xDict:
                    retDictList.append(xDict.dict_value)
                else:
                    retDictList.append(f"{dDesc.strip()}{notFound}")
        return(retDictList)


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
        app_label = 'apputil'
        db_table = 'app_log'

