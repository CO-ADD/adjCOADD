import pandas as pd
import numpy as np
from sequences import Sequence
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
    HEADER_FIELDS = {
        'name':'Name',
        'username':'Username', 
        'first_name':'First Name',  
        'last_name':'Last Name',
        'initials':'Initial',
        'email':'Email',
        'permission':'Permissions',
        'is_appuser':'AppUser',
        'is_active':'Active', 
        }

    username = models.CharField(unique=True, max_length=55, verbose_name='uq user')       # uqjzuegg 
    name = models.CharField(primary_key=True,  max_length=50, verbose_name='user name')          # J.Zuegg
    initials = models.CharField(max_length=5, null=True, blank=True)           # JZG
    organisation = models.CharField(max_length=250, null=True, blank=True)     # University of Queensland
    department = models.CharField(max_length=250, null=True, blank=True)       # Institute for Molecular Bioscience
    group = models.CharField(max_length=50, null=True, blank=True)             # Blaskovich
    phone = models.CharField(max_length=25, null=True, blank=True)             # +61 7 344 62994
    permission = models.CharField(max_length=10, default = 'No', null=False)      # application permissions .. Read, Write, Delete, Admin ..
    is_appuser=models.BooleanField(default=True)

    #------------------------------------------------
    class Meta:
        db_table = 'app_user'
        ordering=['name']

    
    #------------------------------------------------  
    def __str__(self) -> str:
        return f"{self.name}" 

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def get(cls,UserName):
        try:
            retInstance = cls.objects.get(name=UserName)
        except:
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,UserName):
        return cls.objects.filter(name=UserName).exists()

    # --------------------------------------------------------------------------
    def has_permission(self,strPermission) -> bool:
    #
    # Returns True/False if User has strPermission 
    #   checks if strPermission/self.permission in [Read, Write, Delete, Admin]
    # --------------------------------------------------------------------------
        _Permissions = {
            'Read':1,
            'Write':2,
            # 'Delete':3,
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
    # # get field verbose or customized name in the order provided by constants.py

    @classmethod
    def get_fields(cls, fields=None):
        if fields is None:
            fields = cls.HEADER_FIELDS
        if fields:
            select_fields=[fields[f.name] for f in cls._meta.fields if f.name in fields.keys()]
        else:
            select_fields=None
        return select_fields


    # get field name in model Class in the order provided by constants.py
    @classmethod
    def get_modelfields(cls, fields=HEADER_FIELDS):
        if fields:
            model_fields=[f.name for f in cls._meta.fields if f.name in fields.keys()]
        else:
            model_fields=None
        return model_fields

#-------------------------------------------------------------------------------------------------
class AuditModel(models.Model):
    """
    An abstract base class model that provides audit informations 
    """
#-------------------------------------------------------------------------------------------------
    # object status -> indicated by number in database:
    DELETED   = -9
    INVALID   = -1
    UNDEFINED =  0
    VALID     =  1
    CONFIRMED =  2

    OWNER           = "orgdb" # Defaut username
    

    # VALID_STATUS 0: Valid (New or Update) to save, -1: Invalid unable to save, 1: Valid no update required 

    VALID_STATUS    = False
    HEADER_FIELDS   = {}
    CARDS_FIELDS   = {}

    ID_SEQUENCE = "General"
    ID_PREFIX = "AAA"
    ID_PAD = 6

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
        ordering=['-acreated_at']
    
    #------------------------------------------------
    def validate(self,**kwargs):
    #
    # Validates the instance using full_clean
    # 
        retValid = {}
        try:
            self.full_clean(**kwargs)
        except ValidationError as e:
            for key in e.message_dict:
                _field = self._meta.get_field(key)
                #print(f"{key} PK:{_field.primary_key},Null:{_field.null},Blank:{_field.blank}")
                errMsgList = e.message_dict[key]
                retMsg = []
                for errMsg in errMsgList:
                    #print(f"{key} -- {errMsg}")
        
                    if not _field.primary_key:
                        if 'This field cannot be null.' == errMsg: 
                            if not _field.null:
                                retMsg.append(errMsg)
                        elif 'This field cannot be blank.' == errMsg:
                            if not _field.blank:
                                retMsg.append(errMsg)
                        elif 'Ensure that there are no more than' in errMsg:
                            if _field.get_internal_type() != 'DecimalField':
                                retMsg.append(errMsg)
                        else:
                            retMsg.append(errMsg)

                if len(retMsg) > 0 :
                    retValid[key] = "; ".join(retMsg)
                #print(len(e.message_dict[key]))
                # if e.message_dict[key] == ['This field cannot be null.']:
                #     if not self._meta.get_field(key).null:
                #         retValid[key] = ", ".join(e.message_dict[key])
                # elif e.message_dict[key] == ['This field cannot be blank.']:
                #     if not self._meta.get_field(key).blank:
                #         retValid[key] = ", ".join(e.message_dict[key])
                # else:
                #     retValid[key] = ", ".join(e.message_dict[key])
        return(retValid)


    #-------------------------------------------------------------------
    def clean_Fields(self, default_Char="", default_Integer=0, default_Decimal=0.0):
    #
    # Sets 'None' fields in the instance according to Django guidelines 
    #   sets CharField    to "" (empty) or 'default' 
    #   sets IntegerField to 0 or 'default'
    #   sets DecimalField to 0.0 or 'default'
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
            if fType == "DecimalField":
                if hasattr(self,field.name):
                    if getattr(self,field.name) is None:
                        defValue = default_Decimal
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
    def __str__(self) -> str:
        return f"{self.pk}"
    #------------------------------------------------
    def __repr__(self) -> str:
        # return f"{self.__name__}: {self.pk}"
        return f"{self.pk}"

    #------------------------------------------------
    @classmethod
    def get(cls,pkID,verbose=0):
        try:
            retInstance = cls.objects.get(pk=pkID)
        except:
            if verbose:
                print(f"[{cls.__name__} Not Found] {pkID} ")
            retInstance = None
        return(retInstance)
    #------------------------------------------------
    @classmethod
    def exists(cls,pkID,verbose=0):
        return cls.objects.filter(pk=pkID).exists()

   #------------------------------------------------
    @classmethod
    def str_id(cls,clsNo) -> str:
        return(f"{cls.ID_PREFIX}{clsNo:0{cls.ID_PAD}d}")

    #------------------------------------------------
    @classmethod
    def next_id(cls) -> str:
        cls_IDSq=Sequence(cls.ID_SEQUENCE)
        cls_nextNo = next(cls_IDSq)
        cls_strID = cls.str_id(cls_nextNo)
        while cls.objects.filter(pk=cls_strID).exists():
            cls_nextNo = next(cls_IDSq)
            cls_strID = cls.str_id(cls_nextNo)
        return(cls_strID)    

    #------------------------------------------------
    def delete(self,**kwargs):
        appuser=kwargs.get("user")
        kwargs.pop("user",None)
        if appuser is None:
            appuser = ApplicationUser.objects.get(name=self.OWNER)

        self.astatus = self.DELETED
        self.adeleted_id = appuser
        self.adeleted_at = timezone.now()
        super(AuditModel,self).save(**kwargs)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        #
        # Checks for application user
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
    # get field verbose or customized name in the order provided by headerfields
    @classmethod
    def get_fields(cls, fields=None):
        if fields is None:
            fields = cls.HEADER_FIELDS
        if fields:
            fieldsname=[field.name for field in cls._meta.fields]
            select_fields=[fields[f] for f in fields.keys() if f in fieldsname or f.split(".")[0] in fieldsname]
            
        else:
            select_fields=None   
        return select_fields
 
    #------------------------------------------------
    # recurse func to get foreign key related table col value
    def iter_foreignkey(self, nameArray=None, n=None):
        obj_parent=getattr(self, nameArray[0])
        obj = obj_parent
        i=1
        while i < n:
            obj = getattr(obj, nameArray[i])
            i += 1
            if i>=5: # break it !! if it leads to the 5th related table!!
                break
        return obj
    
    #------------------------------------------------
    def get_values(self, fields=None):
        
        from django.db.models import Model
        if fields is None:
            fields = self.HEADER_FIELDS

        value_list=[]
        fieldsname=[field.name for field in self._meta.fields]
        for name in fields.keys():
            n=len(name.split("."))
            nameArray=name.split(".")
            if n>1 and nameArray[0] in fieldsname:
                obj = self.iter_foreignkey(nameArray=nameArray, n=n)            
                if isinstance(fields[name], dict):
                    #if HEADER_FIELDS contains a dictionary for link/url information
                    # for foreignkey link not equal field values
                    url_name = list(list(fields[name].values())[0].keys())[0]
                    if url_name != name:
                        n=len(url_name.split("."))
                        urlArray=url_name.split(".")
                        obj_link= self.iter_foreignkey(nameArray=urlArray, n=n)
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(obj_link))})
                    else:
                        # foreignkey link name equal field values 
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(obj))})
                elif isinstance(obj, Model):
                    value_list.append(obj.pk)   
                elif isinstance(obj, list):
                    varray_to_string=','.join(str(e) for e in obj)
                    value_list.append(varray_to_string)   
                else:
                    value_list.append(obj)
                
            elif name in fieldsname:
                obj=getattr(self, name)
                if obj:
                    # check field value is a dict with link value
                    if isinstance(fields[name], dict):
                        # for slug-name
                        url_name = list(list(fields[name].values())[0].keys())[0]
                        url=getattr(self, url_name)
                        # Append link to the value list
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(url))})
                    elif isinstance(obj, Model):
                        value_list.append(obj.pk)
                    elif isinstance(obj, list):
                        array_to_string=','.join(str(e) for e in obj)
                        value_list.append(array_to_string) 
                    else:   
                        value_list.append(obj)
                else:
                    value_list.append(" ")
                    
        return value_list
        
    # used to get fieds and values for display in Card.html
    def get_fieldsandvalues(self, card_fields=None):
        if card_fields is None:
            card_fields = self.CARDS_FIELDS
        fields=self.__class__.get_fields(fields=card_fields)
        values=self.get_values(fields=card_fields)

        try:
            return [(fields[i], values[i]) for i in range(len(fields))]
        except Exception as err:
            return [("error: ", err)]
    
    # get choice from model field queryset
    @classmethod
    def get_field_choices(cls, field_name = None, filter_dict = None):
        #choice_query = cls.objects.order_by().values_list(field_name).distinct()
        if filter_dict:
            choice_query = cls.objects.filter(**filter_dict).distinct().order_by(field_name).values_list(field_name)            
            #print(f"{field_name} {choice_query.count()} {filter_dict}")
        else:
            choice_query = cls.objects.distinct().order_by(field_name).values_list(field_name)
            #print(f"{field_name} {choice_query.count()} ")
            
        choices = [(i[0], i[0]) for i in choice_query]
        return choices

    



#-------------------------------------------------------------------------------------------------
class Dictionary(AuditModel):
#-------------------------------------------------------------------------------------------------
    HEADER_FIELDS = {
        'dict_value':'Value', 
        'dict_class':'Class',  
        'dict_desc':'Description',
        'dict_sort':'Order',
        'dict_app' :'Application',  
    }
    
    dict_value =models.CharField(primary_key=True, unique=True, max_length=50, verbose_name = "Value"  )
    dict_class= models.CharField(max_length=30, verbose_name = "Class")
    dict_app= models.CharField(max_length=30, verbose_name = "Application")
    dict_desc = models.CharField(max_length=140, blank=True, verbose_name = "Description")
    dict_sort = models.IntegerField(default=0, verbose_name = "Order")
   
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
        return str(self.dict_value)

    def __repr__(self) -> str:
        return f"{self.dict_value}|{self.dict_desc}"
        #return f"[{self.dict_class}] {self.dict_value} ({self.dict_desc})"

    def strtml(self)-> str:
        return f"{self.dict_value} <small class='not-visible'> {self.dict_desc} </small>"

    #------------------------------------------------
    @classmethod
    def get(cls,DictClass,DictValue=None,DictDesc=None,verbose=1):
    #
    # Returns a Dictionary instance if found 
    #    by dict_value
    #    by dict_desc (set dict_value = None)
    #
        if DictValue:
            try:
                retDict = cls.objects.get(dict_value=DictValue, dict_class=DictClass)
            except:
                if verbose:
                    print(f"[Dict Value Not Found] {DictValue} {DictClass}")
                retDict = None
        elif DictDesc:
            try:
                retDict = cls.objects.get(dict_desc=DictDesc, dict_class=DictClass)
            except:
                if verbose:
                    print(f"[Dict Desc Not Found] {DictDesc} {DictClass}")
                retDict = None
        else:
            retDict = None
        return(retDict)

    #------------------------------------------------
    @classmethod
    def exists(cls,DictClass,DictValue=None,DictDesc=None,verbose=1):
    #
    # Returns if Dictionary instance exists
    #    by dict_value
    #    by dict_desc (set dict_value = None)
    #
        if DictValue:
            retValue = cls.objects.filter(dict_value=DictValue, dict_class=DictClass).exists()
        elif DictDesc:
            retValue = cls.objects.filter(dict_desc=DictDesc, dict_class=DictClass).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    @classmethod
    #
    # Returns the objects based on the default filter for dict_class= and aStatus>=0
    #
    def get_filterobj(cls,DictClass,showDeleted=False):
        if showDeleted:
            return cls.objects.filter(dict_class=DictClass)
        else:
            return cls.objects.filter(dict_class=DictClass, astatus__gte=0)

    #------------------------------------------------
    @classmethod
    #
    # Returns Dictionary entries for a DictClass as Choices
    #
    def get_aschoices(cls, DictClass, showDesc = True, sep = " | ", emptyChoice= ('--', ' <empty> ')):
        dictList = None
        choices=(emptyChoice,)
        # comment on initial migrations
        dictList=Dictionary.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc', 'dict_sort')
        if dictList:
            sortedlist = sorted(dictList, key=lambda d: d['dict_sort']) 
            if sortedlist:
                choices_values=tuple([tuple(d.values()) for d in sortedlist])
                if showDesc:
                    choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_values)
                else:
                    choices=tuple((a[0], a[0]) for a in choices_values)
        return choices

    # #------------------------------------------------
    # @classmethod
    # #
    # # Returns Dictionary entries for a DictClass as Choices
    # #
    # def get_aschoices(cls, DictClass, showDesc = True, sep = " | ", emptyChoice= ('--', 'empty')):
    # #def get_Dictionary_asChoice(cls, DictClass, showDesc = True, sep = " | ", emptyChoice= ('--', 'empty')):
    #     dictList=cls.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc', 'dict_sort')
    #     sortedlist = sorted(dictList, key=lambda d: d['dict_sort']) 
    #     if sortedlist:
    #         choices_values=tuple([tuple(d.values()) for d in sortedlist])
    #         if showDesc:
    #             choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_values)
    #         else:
    #             choices=tuple((a[0], a[0]) for a in choices_values)
    #     else:
    #         choices=emptyChoice
    #     return choices
    
    #------------------------------------------------
    #@classmethod
    #
    # Returns Dictionary entries for a DictClass as Choices
    #
    #------------------------------------------------
    @classmethod
    #
    # Returns Dictionary entries for a DictClass as Choices
    #
    def get_DictValues_fromStrList(cls,DictClass,DictValueStr=None,DictDescStr=None,sep=";",notFound="#",verbose=1):
    #-----------------------------------------------------------------------------------
        retDictList = []
        retErrList = []
        if DictValueStr:
            dLst = DictValueStr.split(sep)
            for dVal in dLst:
                xDict = cls.get(DictClass,dVal.strip(),None,verbose=verbose)
                if xDict:
                    retDictList.append(xDict.dict_value)
                else:
                    retErrList.append(dVal.strip())
                    #retDictList.append(f"{dVal.strip()}{notFound}")

        elif DictDescStr:
            dLst = DictDescStr.split(sep)
            for dDesc in dLst:
                xDict = cls.get(DictClass,None,dDesc.strip(),verbose=verbose)
                if xDict:
                    retDictList.append(xDict.dict_value)
                else:
                    retErrList.append(dDesc.strip())
                    #retDictList.append(f"{dDesc.strip()}{notFound}")
        return(retDictList,retErrList)

    # @classmethod
    # def get_DictionaryStrList_asArray(cls,DictClass,DictValueStr=None,DictDescStr=None,sep=";",notFound="#"):
    # #-----------------------------------------------------------------------------------
    #     retDictList = []
    #     if DictValueStr:
    #         dLst = DictValueStr.split(sep)
    #         for dVal in dLst:
    #             xDict = cls.get(DictClass,dVal.strip(),None)
    #             if xDict:
    #                 retDictList.append(xDict.dict_value)
    #             else:
    #                 retDictList.append(f"{dVal.strip()}{notFound}")
    #     elif DictDescStr:
    #         dLst = DictDescStr.split(sep)
    #         for dDesc in dLst:
    #             xDict = cls.get(DictClass,None,dDesc.strip())
    #             if xDict:
    #                 retDictList.append(xDict.dict_value)
    #             else:
    #                 retDictList.append(f"{dDesc.strip()}{notFound}")
    #     return(retDictList)


#-------------------------------------------------------------------------------------------------
class ApplicationLog(models.Model):
#-------------------------------------------------------------------------------------------------
    HEADER_FIELDS = {
        'log_code':'Code', 
        'log_proc':'Procedure',  
        'log_type':'Type',
        'log_time':'Time',
        'log_status':'Status',
        'log_user':'User',
        'log_object' :'Object',  
        'log_desc' :'Description',  
    }

    OWNER     = "orgdb"

    log_code = models.CharField(max_length=15, blank=True, db_index = True, editable=False,verbose_name = "Code")
    log_proc = models.CharField(max_length=50, blank=True, db_index = True, editable=False,verbose_name = "Procedure")
    log_type = models.CharField(max_length=15, blank=True, db_index = True, editable=False,verbose_name = "Type")
    log_time = models.DateTimeField(auto_now=True, editable=False,verbose_name = "Time")
    log_user = models.ForeignKey(ApplicationUser, blank=True, verbose_name = "User", on_delete=models.DO_NOTHING, 
        db_column="log_user", related_name="%(class)s_user")
    log_object = models.CharField(max_length=150, blank=True, db_index = True, editable=False,verbose_name = "Object")
    log_desc = models.CharField(max_length=1024, blank=True, editable=False,verbose_name = "Description")
    log_status = models.CharField(max_length=15, blank=True, db_index = True, editable=False,verbose_name = "Status")

    class Meta:
        app_label = 'apputil'
        db_table = 'app_log'
        ordering=['log_time','log_type']
        indexes = [
            models.Index(name="log_code_idx",fields=['log_code']),
            models.Index(name="log_proc_idx",fields=['log_proc']),
            models.Index(name="log_type_idx",fields=['log_type']),
            models.Index(name="log_status_idx",fields=['log_status']),
            models.Index(name="log_object_idx",fields=['log_object']),
        ]

    #------------------------------------------------
    @classmethod
    #
    # Saves an Log Entry
    #
    def add(cls, LogCode, LogProc, LogType,  LogUser, LogObject, LogDesc, LogStatus):
        log_inst = cls()
        log_inst.log_code = LogCode
        log_inst.log_proc = LogProc
        log_inst.log_type = LogType
        if LogUser is None:
            LogUser = ApplicationUser.get(cls.OWNER)
        log_inst.log_user = LogUser
        log_inst.log_object = LogObject
        log_inst.log_desc = LogDesc
        log_inst.log_status = LogStatus
        log_inst.save()

    #------------------------------------------------
    # get field verbose or customized name in the order provided by headerfields
    @classmethod
    def get_fields(cls, fields=None):
        if fields is None:
            fields = cls.HEADER_FIELDS
        if fields:
            fieldsname=[field.name for field in cls._meta.fields]
            select_fields=[fields[f] for f in fields.keys() if f in fieldsname or f.split(".")[0] in fieldsname]
            
        else:
            select_fields=None   
        return select_fields

    #------------------------------------------------
    # recurse func to get foreign key related table col value
    def iter_foreignkey(self, nameArray=None, n=None):
        obj_parent=getattr(self, nameArray[0])
        obj = obj_parent
        i=1
        while i < n:
            obj = getattr(obj, nameArray[i])
            i += 1
            if i>=5: # break it !! if it leads to the 5th related table!!
                break
        return obj
    
    #------------------------------------------------
    def get_values(self, fields=None):
        
        from django.db.models import Model
        if fields is None:
            fields = self.HEADER_FIELDS

        value_list=[]
        fieldsname=[field.name for field in self._meta.fields]
        for name in fields.keys():
            n=len(name.split("."))
            nameArray=name.split(".")
            if n>1 and nameArray[0] in fieldsname:
                obj = self.iter_foreignkey(nameArray=nameArray, n=n)            
                if isinstance(fields[name], dict):
                    #if HEADER_FIELDS contains a dictionary for link/url information
                    # for foreignkey link not equal field values
                    url_name = list(list(fields[name].values())[0].keys())[0]
                    if url_name != name:
                        n=len(url_name.split("."))
                        urlArray=url_name.split(".")
                        obj_link= self.iter_foreignkey(nameArray=urlArray, n=n)
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(obj_link))})
                    else:
                        # foreignkey link name equal field values 
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(obj))})
                elif isinstance(obj, Model):
                    value_list.append(obj.pk)   
                elif isinstance(obj, list):
                    varray_to_string=','.join(str(e) for e in obj)
                    value_list.append(varray_to_string)   
                else:
                    value_list.append(obj)
                
            elif name in fieldsname:
                obj=getattr(self, name)
                if obj:
                    # check field value is a dict with link value
                    if isinstance(fields[name], dict):
                        # for slug-name
                        url_name = list(list(fields[name].values())[0].keys())[0]
                        url=getattr(self, url_name)
                        # Append link to the value list
                        value_list.append({obj: list(list(fields[name].values())[0].values())[0].replace('{VALUE1}', str(url))})
                    elif isinstance(obj, Model):
                        value_list.append(obj.pk)
                    elif isinstance(obj, list):
                        array_to_string=','.join(str(e) for e in obj)
                        value_list.append(array_to_string) 
                    else:   
                        value_list.append(obj)
                else:
                    value_list.append(" ")
                    
        return value_list
        
    # get choice from model field queryset
    @classmethod
    def get_field_choices(cls, field_name = None):
        choice_query = cls.objects.order_by().values_list(field_name).distinct()
        choices = [(i[0], i[0]) for i in choice_query]
        return choices

# #-------------------------------------------------------------------------------------------------
# class Image(AuditModel):
# #-------------------------------------------------------------------------------------------------
#     HEADER_FIELDS = {
#         'image_name':'Name', 
#         'image_file':'Image',  
#         'image_type':'Type',  
#         'image_desc':'Description',
#         'image_source':'Source',
#         'image_object':'Object',
#         'image_objectid':'Object ID',
#     }
    
#     image_name =models.CharField(max_length=120,  unique=True, verbose_name = "Name")
#     image_file= models.ImageField(upload_to='images/', verbose_name = "Image")
#     image_type = models.CharField(max_length=25, verbose_name = "Type")
#     image_desc = models.CharField(max_length=140, blank=True, verbose_name = "Description", default = "Description")
#     image_source = models.CharField(max_length=50, blank=True, verbose_name = "Source", default = "source")

#     class Meta:
#         app_label = 'apputil'
#         db_table = 'app_image'
#         ordering=['image_name',]
#         indexes = [
#             models.Index(name="img_name_idx",fields=['image_name']),
#             models.Index(name="img_scr_idx",fields=['image_source']),
#         ]

#     #------------------------------------------------
#     def __str__(self) -> str:
#         return str(self.image_name)

#     def __repr__(self) -> str:
#         return f"[{self.image_name}] {self.image_object} ({self.image_objectid})"

#     #------------------------------------------------
#     @classmethod
#     def get(cls,ImgName,ImgObj=None,ImgObjID=None,verbose=1):
#     #
#     # Returns a Image instance if found 
#     #    by ImgName
#     #    by ImgObj & ImgObjID
#     #
#         if ImgName:
#             try:
#                 retDict = cls.objects.get(image_name=ImgName)
#             except:
#                 if verbose:
#                     print(f"[Image Not Found] {ImgName}")
#                 retDict = None
#         elif ImgObj:
#             try:
#                 retDict = cls.objects.get(image_object=ImgObj, image_objectid=ImgObjID)
#             except:
#                 if verbose:
#                     print(f"[Image Not Found] {ImgObj} {ImgObjID}")
#                 retDict = None
#         else:
#             retDict = None
#         return(retDict)

#     #------------------------------------------------
#     @classmethod
#     def exists(cls,ImgName,ImgObj=None,ImgObjID=None,verbose=1):
#     #
#     # Returns if Image instance exists
#     #    by ImgName
#     #    by ImgObj & ImgObjID
#     #
#         if ImgName:
#             retValue = cls.objects.filter(image_name=ImgName).exists()
#         elif ImgObj:
#             retValue = cls.objects.filter(image_object=ImgObj, image_objectid=ImgObjID).exists()
#         else:
#             retValue = False
#         return(retValue)

#-------------------------------------------------------------------------------------------------
class Document(AuditModel):
#-------------------------------------------------------------------------------------------------
    HEADER_FIELDS = {
        'doc_name':'Name', 
        'doc_file':'Document',  
        'doc_type':'Type',  
        'doc_desc':'Description',
        'doc_source':'Source',
    }
    ID_SEQUENCE = 'Document'
    ID_PREFIC = 'DOC'
    ID_PAD = 5
    
    doc_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Doc ID")
    doc_name =models.CharField(max_length=120, blank=False, verbose_name = "Name"  )
    doc_file= models.FileField(upload_to='documents/', verbose_name = "Document")
    doc_type = models.CharField(max_length=25, verbose_name = "Type")
    doc_desc = models.CharField(max_length=140, blank=True, verbose_name = "Description", default = "description")
    doc_source = models.CharField(max_length=50, blank=True, verbose_name = "Source", default = "source")

    class Meta:
        app_label = 'apputil'
        db_table = 'app_document'
        ordering=['doc_name',]
        indexes = [
            models.Index(name="doc_name_idx",fields=['doc_name']),
            models.Index(name="doc_scr_idx",fields=['doc_source']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"[{self.doc_id}] {self.doc_name} ({self.doc_file})"

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.doc_id:
            self.doc_id = self.next_id()
            if self.doc_id: 
                super(Document, self).save(*args, **kwargs)
        else:
            super(Document, self).save(*args, **kwargs) 
