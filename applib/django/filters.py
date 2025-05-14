"""
Used and import for all application filter views
"""
from datetime import datetime
import pandas as pd

from django import forms
from django.views.generic import ListView
from django.db.models import Q, CharField, TextField, ForeignKey, IntegerField, Func, Value, ManyToManyField
from django.db.models import CharField, Value
from django.db.models.functions import Cast
from django.contrib.postgres.fields import ArrayField
from django.db.models.expressions import RawSQL
from django.contrib.postgres.search import TrigramSimilarity
from django.db.models.functions import Greatest
from django.core.validators import MinLengthValidator

from django_filters import FilterSet, CharFilter, ChoiceFilter

# -- create a function for search all fields--
#--------------------------------------------------------------------------
def get_all_fields_q_object(model, search_value, exclude_fields=None, prefix=None, submodel=False):
#--------------------------------------------------------------------------
    
    q_object = Q()
    exclude_fields = exclude_fields or []
    searchfields=[]
    if submodel:
        searchfields.append(model._meta.pk)
    else:
        searchfields=[field for field in model._meta.get_fields()]
    for field in searchfields:
        if field.name in exclude_fields:
            continue
        lookup_field_name = f"{prefix}__{field.name}" if prefix else field.name   
        if isinstance(field, (CharField, TextField)):
            q_object |= Q(**{f"{lookup_field_name}__icontains": search_value})        
        elif isinstance(field, ForeignKey):
            related_model = field.related_model
            related_q_object = get_all_fields_q_object(related_model, search_value, exclude_fields=exclude_fields, prefix=lookup_field_name, submodel=True)
            q_object |= related_q_object
        elif isinstance(field, IntegerField):
            try:
                int_value = int(search_value)
                q_object |= Q(**{lookup_field_name: int_value})
            except ValueError:
                pass       
        elif isinstance(field, ArrayField):
            q_object |= Q(**{f'{lookup_field_name}__icontains': search_value})       
        elif isinstance(field, ManyToManyField):
            related_model = field.related_model
            related_q_object = get_all_fields_q_object(related_model, search_value, exclude_fields=exclude_fields, prefix=lookup_field_name, submodel=True)
            q_object |= related_q_object
        # # Add more field types as needed...

    return q_object

#--------------------------------------------------------------------------
def get_all_fields_q_object_deep(model, search_value, exclude_fields=None, prefix=None):
#--------------------------------------------------------------------------

    q_object = Q()
    exclude_fields = exclude_fields or []

    for field in model._meta.get_fields():
        if field.name in exclude_fields:
            continue
        lookup_field_name = f"{prefix}__{field.name}" if prefix else field.name
        if isinstance(field, (CharField, TextField)):
            q_object |= Q(**{f"{lookup_field_name}__icontains": search_value})
        
        elif isinstance(field, ForeignKey):
            related_model = field.related_model
            related_q_object = get_all_fields_q_object_deep(related_model, search_value, exclude_fields=exclude_fields, prefix=lookup_field_name)
            q_object |= related_q_object
        elif isinstance(field, IntegerField):
            try:
                int_value = int(search_value)
                q_object |= Q(**{lookup_field_name: int_value})
            except ValueError:
                pass       
        elif isinstance(field, ArrayField):
            q_object |= Q(**{f'{lookup_field_name}__icontains': search_value})        
        elif isinstance(field, ManyToManyField):
            related_model = field.related_model
            related_q_object = get_all_fields_q_object_deep(related_model, search_value, exclude_fields=exclude_fields, prefix=lookup_field_name)
            q_object |= related_q_object
        # # Add more field types as needed...

    return q_object

# -- Filterset base Class--
# from adjcoadd.constants import CharToChoice_filterList
# from django.contrib import messages

#--------------------------------------------------------------------------
class Base_Filter(FilterSet):
#--------------------------------------------------------------------------

    Search_all_fields = CharFilter(method='filter_all_fields', 
                                widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields', 'minlength':'3' }), 
                                validators=[MinLengthValidator(3)])

    def __init__(self, *args, **kwargs):
        deep=kwargs.pop('deep') # switcher deep search or one-table search
        filterset_dict=kwargs.pop('filterset_dict',None) # switcher deep search or one-table search
        
        super().__init__(*args, **kwargs)

        if deep == True:
            self.filters['Search_all_fields'].method = 'filter_all_fields_deep'
        for field in self.filters:
            if 'CharFilter' == self.filters[field].__class__.__name__:
                self.filters[field].lookup_expr='icontains'
        
        # Loop through all fields, find Charfield to Choicefield
        # for field in self.filters:
        #     if str(field) in CharToChoice_filterList:
        #         self.filters[str(field)].extra["choices"] = ChoiceFilter(choices=self.Meta.model.get_field_choices(field_name=str(field)))

    def update_choice_filters(self):
        print(f"Filterbase_base.update_choice_filters")
    
    def multichoices_filter(self, queryset, name, value):
        lookup='__'.join([name, 'overlap'])
        return queryset.filter(**{lookup: value})
    
    def filter_all_fields(self, queryset, name, value):
        if value:
            exclude_fields = ['password','astatus',]
            q_object = get_all_fields_q_object(self._meta.model, value, exclude_fields=exclude_fields)
            return queryset.filter(q_object)
        return queryset
    
    def filter_all_fields_deep(self, queryset, name, value):
        if value:
            exclude_fields = ['password','astatus',]
            q_object = get_all_fields_q_object_deep(self._meta.model, value, exclude_fields=exclude_fields)
            return queryset.filter(q_object)
        return queryset
    
         
#--------------------------------------------------------------------------
class BaseStatus_Filter(Base_Filter):
#--------------------------------------------------------------------------

    # Filter primary queryset for valid (not deleted) entries with aStatus >= 0 
    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)


    
# # utils for filteredListView method def ordered_by
# #--------------------------------------------------------------------------
# def find_item_index(lst, item):
# #--------------------------------------------------------------------------

#     for i, element in enumerate(lst):
#         if isinstance(element, dict):
#             if item in element.keys():
#                 return i
#         elif element == item:
#             return i
        
#     return -1

# --Filter view base class--

    
    # 
 