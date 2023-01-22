import os
import re
import unicodedata
import django_filters
import pandas as pd

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction
from django.views.generic import ListView

from .models import Dictionary
from adjcoadd.constants import *



# ===================================Dictionary query convert to choice Tuples========================================================================#

def get_DictonaryChoices_byDictClass(ModelName, DictClass, sep='|'):
    options=ModelName.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc')
    if options:
        choices_values=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_values)
    else:
        choices=(('--', 'empty'),)
    return choices
# ------------------------Only use dict_value
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
# Base Class for all models list/card view
class FilteredListView(ListView):
    filterset_class = None
    paginate_by=50
    model_fields=None
    order_by=None

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
        print("no order")
        return self.filterset.qs.distinct()

    def get_context_data(self,**kwargs):
        context = super().get_context_data(**kwargs)
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

    

    