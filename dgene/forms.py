import django_filters
from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField

from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import Filterbase
from .models import Gene

class Gene_form(ModelForm):

     
    gene_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class="gene_type"), required=False)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False,)
   
    def __init__(self, *args, **kwargs):    
        super().__init__(*args, **kwargs)
        self.fields['gene_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Gene.Choice_Dictionary['gene_type'], astatus__gte=0)]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs

    def create_field_groups(self):
        self.group1 = [self[name] for name in Gene.FORM_GROUPS['Group_gene']]
        self.group2 = [self[name] for name in Gene.FORM_GROUPS['Group_protein']]
 
    
    class Meta:
        model=Gene
        exclude = ['gene_id']
 
# -------------fitlerset Forms---------------------------------------------------------------

class Genefilter(Filterbase):
   
    Gene_Type=django_filters.ChoiceFilter(field_name='gene_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Gene_name=django_filters.CharFilter(field_name='gene_name', lookup_expr='icontains', label='Gene Name')
    Othername=django_filters.CharFilter(field_name='gene_othernames', lookup_expr='icontains', label='Other Name')
    Gene_Class=django_filters.CharFilter(field_name='protein_class', lookup_expr='icontains', label='Gene Class')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Gene_Type"].extra['choices']=[(obj.dict_value, obj) for obj in Dictionary.objects.filter(dict_class=Gene.Choice_Dictionary['gene_type'], astatus__gte=0)]
        
    class Meta:
        model=Gene
        fields=[ 'Gene_Type']
