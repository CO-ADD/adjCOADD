from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter
from dorganism.models import  Taxonomy

from apputil.models import Dictionary, ApplicationUser, Document
from applib.django.base.filters import BaseStatus_Filter
from adjcoadd.constants import ORGANISM_CLASSES, CELL_CLASSES 
from dpeptide.models import Peptide, Peptide_Batch

#=================================================================================================
# Peptide
#=================================================================================================
class Peptide_Filter(BaseStatus_Filter):
    
    peptide_id = CharFilter(field_name='peptide_id', lookup_expr='icontains')
    peptide_name = CharFilter(field_name='peptide_name', lookup_expr='icontains')
    peptide_notes = CharFilter(field_name='peptide_notes', lookup_expr='icontains')
    peptide_type = MultipleChoiceFilter(field_name='peptide_type', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    peptide_panel = MultipleChoiceFilter(field_name='peptide_panel', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])

    mta_status = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Peptide.DICTIONARY_FIELDS['mta_status'], astatus__gte=0))
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["peptide_type"].extra["choices"]=Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_type'], showDesc = False)
        self.filters["peptide_panel"].extra["choices"]=Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_panel'], showDesc = False)

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=Peptide
        fields=[ 'peptide_id', 'peptide_name','peptide_notes', 'peptide_type', 'peptide_panel', 'mta_status' ]

# -----------------------------------------------------------------
class Peptide_CreateForm(forms.ModelForm):

    peptide_name= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    peptide_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    # prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    #organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    #collect_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
   
    def __init__(self, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        super(Peptide_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['peptide_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_type'], showDesc=False),)
        self.fields['peptide_type'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true',})

        self.fields['peptide_panel'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Peptide.DICTIONARY_FIELDS['peptide_panel'], showDesc=False),)
        self.fields['peptide_panel'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true'})

        self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Peptide.DICTIONARY_FIELDS['mta_status'])]
        self.create_field_groups()

        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    # def clean_organism_name(self):       
    #     print(self.organism_name)
    #     data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
    
    #     if data:
    #         if str(data.org_class) in CELL_CLASSES:
    #             return data
    #         else:
    #             self.add_error('org_name', 'Create failed with invalid Organism class')
                
    #     else:
    #         self.add_error('org_name', "Found No Organism")

    #     return data            


    def create_field_groups(self):
        if len(Peptide.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Peptide.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   
    
    class Meta:
        model=Peptide
        exclude=['peptide_id',  'assoc_documents'] 

# -----------------------------------------------------------------
class Peptide_UpdateForm(Peptide_CreateForm):     
    class Meta:
        model=Peptide
        exclude=['peptide_id', 'assoc_documents'] 

