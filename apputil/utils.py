import os
from pathlib import Path
from datetime import datetime
import re
import unicodedata
import django_filters
import pandas as pd
from asgiref.sync import sync_to_async
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction
from django.views.generic import ListView
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.conf import settings

from .models import Dictionary
from adjcoadd.constants import *

# ==========utilized in Create User folder under /uploads ==========
def file_location(req):
    location=settings.MEDIA_ROOT+'/'+str(req.user)
    return location

# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")

# ==========Super UserRequire Mixin===================================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = '/'

    def test_func(self):
        return self.request.user.has_permission('Admin')

# ============Override filename in FileStorage========================================
from django.core.files.storage import FileSystemStorage

class OverwriteStorage(FileSystemStorage):
    
    def get_available_name(self, name, max_length=None):
        """
        Returns a filename that's free on the target storage system, and
        available for new content to be written to.
        """
        # If the filename already exists, remove it as if it was a true file system
    
        if self.exists(name):
            os.remove(os.path.join(self.location, name))
        return name

# Model Validation utilities  =====================================================================================

class Validation_Log():

    def __init__(self,logProcess,logTypes= ['Error','Warning','Info']):
        self.logProcess = logProcess
        self.logTypes = logTypes
        self.nLogs = {}
        self.Logs={}

        for t in self.logTypes:
            self.nLogs[t] = 0
            self.Logs[t] = []
           

    def add_log(self, logType, logDesc, logItem, logHelp):
        lDict = {
            'Process': self.logProcess, 
            'Description': logDesc, 
            'Item': str(logItem), 
            'Help': logHelp,
            'Time': datetime.now() }
        logType = logType[0].upper()+logType[1:].lower()
        if logType in self.logTypes:
            self.Logs[logType].append(lDict)
            self.nLogs[logType] = self.nLogs[logType] + 1
        
    def show(self,logTypes= ['Error','Warning', 'Info']):
        info={} #info=[]
        for t in logTypes:
            # print(f"-- {t.upper():8} ({self.nLogs[t]:3}) ------------------------------------------------------")
            info[t]=[]
            for l in self.Logs[t]:
                print(f"{l['Process']}-{l['Description']} ({l['Item']}) {l['Help']} ")
                description=str(l['Description']).replace("'", "").replace('"', '')
                print_info=f"{l['Process']}_{description}_{l['Item']}_{l['Help']}"
                info[t].append(print_info) # info.append(print_info)
       
        return info



    
#-----------------------------------------------------------------------------------
def instance_dict(instance, key_format=None):
    "Returns a dictionary containing field names and values for the given instance"
#-----------------------------------------------------------------------------------
    from django.forms.models import model_to_dict
    model_to_dict(instance, fields=[field.name for field in instance._meta.fields]) 

# ===================================Dictionary query convert to choice Tuples========================================================================#

def get_DictonaryChoices_byDictClass(ModelName, DictClass, sep='|'):
    options=ModelName.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc')
    if options:
        choices_values=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_values)
    else:
        choices=(('--', 'empty'),)
    return choices
    
# ------------------------Only use dict_value----------------
def get_DictonaryChoicesValue_byDictClass(ModelName, DictClass, sep='|'):
    options=ModelName.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc')
    if options:
        choices_values=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]) for a in choices_values)
    else:
        choices=(('--', 'empty'),)
    return choices

#-----------------------------------------------------------------------------------
def slugify(value, lower=False, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize("NFKC", value)
    else:
        value = (
            unicodedata.normalize("NFKD", value)
            .encode("ascii", "ignore")
            .decode("ascii")
        )
    value = re.sub(r"[^\w\s-]", "", value)
    value = re.sub(r"[-\s]+", "-", value).strip("-_")

    if lower:
        return value.lower()
    else:
        return value
#-----------------------------------------------------------------------------------
 


#  #####################Django Filter View#################
class Filterbase(django_filters.FilterSet):
   
    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)
   
    def multichoices_filter(self, queryset, name, value):
        lookup='__'.join([name, 'overlap'])
        return queryset.filter(**{lookup: value})




# Base Class for all models list/card view
class FilteredListView(ListView):
    filterset_class = None
    paginate_by=50
    model_fields=None
    order_by=None
    context_list=''

   
    def get_queryset(self):
        # Get the queryset however you usually would.  For example:
        queryset = super().get_queryset()
        # Then use the query parameters and the queryset to
        # instantiate a filterset and save it as an attribute
        # on the view instance for later.
        self.filterset = self.filterset_class(self.request.GET, queryset=queryset)
        # Return the filtered queryset
        order=self.get_order_by()
        if order:
            return self.filterset.qs.distinct().order_by(order)
        # print(self.filterset.qs.distinct())
        return self.filterset.qs.distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        self.context_list=context['object_list']
        # print(context)
        # Pass the filterset to the template - it provides the form.
        context['filter'] = self.filterset
        context['paginate_by']=self.get_paginate_by(self, **kwargs)
        context['fields']=self.model.get_fields(fields=self.model_fields)
        context['model_fields']=self.model.get_modelfields(fields=self.model_fields)
      
        return context

    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)
        return paginate_by
 
    def get_order_by(self):
       
        order_by=self.request.GET.get("order_by", self.order_by) or None
        model_constants_field=self.model_fields
        acs_decs=""
        if order_by:
            order_field=""
            if order_by[0]=="-":
                acs_decs=order_by[0]
                order_field=order_by[1:]
            else:
                order_field=order_by
                
            if order_field in model_constants_field.values():
                order_by=acs_decs+ list(model_constants_field.keys())[list(model_constants_field.values()).index(order_field)]
           
            return order_by
        return order_by


a=[tuple(d.values()) for d in Dictionary.objects.order_by().values('dict_class').distinct()]
choice_class=[(x[0], x[0]) for x in a]
class Dictionaryfilter(Filterbase):
    dict_class = django_filters.ChoiceFilter(choices=choice_class)
    #   dict_value = django_filters.CharFilter(lookup_expr='icontains')
    #   dict_desc = django_filters.CharFilter(lookup_expr='icontains')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['dict_class'].label='Class'
      
    class Meta:
        model=Dictionary
        fields=['dict_class']


#file path on server:

# Define full path name
def get_filewithpath( file_name=None):
    if settings.DEVELOPMENT:
        file_path = f"static/images/{file_name}.svg"
   
    else:
        Base_dir = Path(__file__).resolve().parent.parent.parent
        FILES_DIR=os.path.abspath(os.path.join(Base_dir, 'static/images'))
        file_path=os.path.join(FILES_DIR, f"{file_name}.svg")
    return file_path