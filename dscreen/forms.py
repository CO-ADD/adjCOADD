from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter

from apputil.models import Dictionary, ApplicationUser, Document
from apputil.utils.filters_base import Filterbase
#from adjcoadd.constants import ORGANISM_CLASSES


#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50

#=================================================================================================
# Screen_Run
#=================================================================================================
class ScreenRun_Filter(Filterbase):
    
    # run_id = CharFilter(field_name='run_id', lookup_expr='icontains')
    # run_name = CharFilter(field_name='run_name', lookup_expr='icontains')
    # run_type = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    # run_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    # Strain = CharFilter(field_name='strain_ids', lookup_expr='icontains')
    # Notes = CharFilter(field_name='strain_notes', lookup_expr='icontains')
    # Type = MultipleChoiceFilter(field_name='run_type', method='multichoices_filter', 
    #                                          widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    # Panel = MultipleChoiceFilter(field_name='strain_panel', method='multichoices_filter', 
    #                                          widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    # MTA = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Organism.DICTIONARY_FIELDS['mta_status'], astatus__gte=0))
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.fields['run_type'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        # self.fields['run_status'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]
        # self.filters["Type"].extra["choices"]=Dictionary.get_aschoices(Organism.DICTIONARY_FIELDS['strain_type'], showDesc = False)
        # self.filters["Panel"].extra["choices"]=Dictionary.get_aschoices(Organism.DICTIONARY_FIELDS['strain_panel'], showDesc = False)
        # for i in self.filters:
        #     self.filters[i].label=i
   
    class Meta:
        model=Screen_Run
        fields=[ 'run_id', 'run_name','run_type','run_status']

# -----------------------------------------------------------------
class ScreenRun_CreateForm(forms.ModelForm):

    def __init__(self, organism_name=None, *args, **kwargs): 
        super(ScreenRun_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['run_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Screen_Run.DICTIONARY_FIELDS['run_type'], showDesc=False),)
        self.fields['run_status'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Screen_Run.DICTIONARY_FIELDS['run_status'], showDesc=False),)

    class Meta:
        model=Screen_Run
        exclude=Screen_Run.CALCULATED_FIELDS 

# -----------------------------------------------------------------
class ScreenRun_UpdateForm(ScreenRun_CreateForm):     
    class Meta:
        model=Screen_Run
        exclude=Screen_Run.CALCULATED_FIELDS  
