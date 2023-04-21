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
from dorganism.models import Organism_Batch
from .models import Gene, ID_Pub, ID_Sequence, WGS_FastQC, WGS_CheckM

# --Gene Forms--
## 
class Gene_form(ModelForm):
    gene_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class="gene_type"), required=False)
   
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
 
## fitler forms
class Genefilter(Filterbase):
   
    gene_type=django_filters.ChoiceFilter(choices=[])
    # Gene_name=django_filters.CharFilter(field_name='gene_name', lookup_expr='icontains', label='Gene Name')
    # Othername=django_filters.CharFilter(field_name='gene_othernames', lookup_expr='icontains', label='Other Name')
    # Gene_Class=django_filters.CharFilter(field_name='protein_class', lookup_expr='icontains', label='Gene Class')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["gene_type"].extra['choices']=[(obj.dict_value, obj) for obj in Dictionary.objects.filter(dict_class=Gene.Choice_Dictionary['gene_type'], astatus__gte=0)]
        
    class Meta:
        model=Gene
        fields=list(model.HEADER_FIELDS.keys())

# --ID_Pub Forms--
## 
class ID_Pub_form(ModelForm):
    id_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class="id_type"), required=False)
   
    def __init__(self, *args, **kwargs):    
        super().__init__(*args, **kwargs)
        self.fields['id_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=ID_Pub.Choice_Dictionary['id_type'], astatus__gte=0)]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs

    def create_field_groups(self):
        pass
 
    
    class Meta:
        model=ID_Pub
        fields="__all__"
 
## fitler forms
class ID_Pubfilter(Filterbase):
   
    # id_type=django_filters.ChoiceFilter(choices=[])
    # Gene_name=django_filters.CharFilter(field_name='gene_name', lookup_expr='icontains', label='Gene Name')
    # Othername=django_filters.CharFilter(field_name='gene_othernames', lookup_expr='icontains', label='Other Name')
    # Gene_Class=django_filters.CharFilter(field_name='protein_class', lookup_expr='icontains', label='Gene Class')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.filters["id_type"].extra['choices']=[(obj.dict_value, obj) for obj in Dictionary.objects.filter(dict_class=ID_Pub.Choice_Dictionary['id_type'], astatus__gte=0)]
        
    class Meta:
        model=ID_Pub
        fields=list(model.HEADER_FIELDS.keys())

# --ID_Sequence forms--
##
class Sequence_form(ModelForm):

    id_organisms=SimpleArrayField(forms.CharField(), required=False)
    id_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class="id_type"), required=False)
    orgbatch_id=forms.ModelChoiceField(queryset=Organism_Batch.objects.all(), required=False,)
    id_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
   
    def __init__(self, *args, **kwargs):    
        super().__init__(*args, **kwargs)
        self.fields['id_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=ID_Sequence.Choice_Dictionary['id_type'], astatus__gte=0)]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs

    def create_field_groups(self):
        pass
  
    class Meta:
        model=ID_Sequence
        fields="__all__"

## fitler forms
class Sequencefilter(Filterbase):
    run_id=django_filters.CharFilter(lookup_expr='icontains', label='RUN ID')
    id_organisms=django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=[] )
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    class Meta:
        model=ID_Sequence
        fields=list(model.HEADER_FIELDS.keys())

# WGS_FastQCfilter forms
class FastQCfilter(Filterbase):
    Run_Id=django_filters.CharFilter(field_name='run_id', lookup_expr='icontains', label='RUN ID')
    # OrgBatch_ID=django_filters.ModelChoiceFilter(field_name='orgbatch_id', queryset=Organism_Batch.objects.filter(astatus__gte=0))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model=WGS_FastQC
        fields=list(model.HEADER_FIELDS.keys())#[ 'Run_Id']

# WGS_CheckMfilter forms
class CheckMfilter(Filterbase):
    Run_Id=django_filters.CharFilter(field_name='run_id', lookup_expr='icontains', label='RUN ID')
    # OrgBatch_ID=django_filters.ModelChoiceFilter(field_name='orgbatch_id', queryset=Organism_Batch.objects.filter(astatus__gte=0))
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.create_field_groups()
        # self.group1 = [self.filters[name] for name in list(WGS_CheckM.HEADER_FIELDS.keys())]
    
    def create_field_groups(self):
        self.group1 = [self.filters[name] for name in list(WGS_CheckM.HEADER_FIELDS.keys())]
        print(self.group1[0].label)

    class Meta:
        model=WGS_CheckM
        fields=list(model.HEADER_FIELDS.keys()) #[ 'Run_Id']