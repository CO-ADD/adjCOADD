#import django_filters
from django import forms
from django.db import models
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import CharFilter, ChoiceFilter

from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import Filterbase
from dorganism.models import Organism_Batch
from dgene.models import Genome_Sequence, ID_Pub, ID_Sequence, WGS_FastQC, WGS_CheckM, Gene, AMR_Genotype  


#=================================================================================================
# Genome Sequences
#=================================================================================================
class GenomeSeq_Filter(Filterbase):

    f_OrgBatchID = CharFilter(field_name='orgbatch_id__orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model=Genome_Sequence
        fields = ['f_OrgBatchID','f_OrgName']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.orgbatch_id',
                   'orgbatch_id.organism_id.organism_name',
                   ]

class GenomeSeq_Form(ModelForm):
    def __init__(self, *args, **kwargs):    
        super().__init__(*args, **kwargs)

    class Meta:
        model=Genome_Sequence
        fields= ['seq_id']

#=================================================================================================
# ID_Pub - Identification from Public sources
#=================================================================================================
# --ID_Pub Forms--
## 
class IDPub_Form(ModelForm):
    #id_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class="id_type"), required=False)
   
    def __init__(self, *args, **kwargs):    
        super().__init__(*args, **kwargs)
        # self.fields['id_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=ID_Pub.Choice_Dictionary['id_type'], astatus__gte=0)]
        # self.create_field_groups()
        # for field in self.fields.values():
        #     if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
        #         # Add the 'group-input' class to the widget attrs
        #         attrs = field.widget.attrs
        #         attrs['class'] = attrs.get('class', '') + 'input-group'
        #         field.widget.attrs = attrs

    # def create_field_groups(self):
    #     pass
     
    class Meta:
        model=ID_Pub
        fields= ['id_type']
 
## fitler forms
class IDPub_Filter(Filterbase):
    #id_organisms=CharFilter(method='filter_arrayfields')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model=ID_Pub
        fields= ['id_type']
        #fields=list(model.HEADER_FIELDS.keys())

#=================================================================================================
# ID_Seq - Identification from Sequence
#=================================================================================================
class IDSeq_Filter(Filterbase):
    f_OrgBatchID = CharFilter(field_name='orgbatch_id__orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
    kraken_organisms = CharFilter(field_name='kraken_organisms', lookup_expr='icontains',label="Kraken2 Organisms")

    #id_organisms=MultipleChoiceFilter(method='multichoices_filter', choices=[] )
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model = ID_Sequence
        fields = ['f_OrgBatchID','f_OrgName']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.orgbatch_id',
                   'orgbatch_id.organism_id.organism_name',
                   ]

class IDSeq_Form(ModelForm):

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

#=================================================================================================
# WGS_FastQC - FastQ QC
#=================================================================================================
class WGS_FastQC_Filter(Filterbase):
    f_OrgBatchID = CharFilter(field_name='orgbatch_id__orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class Meta:
        model=WGS_FastQC
        fields = ['f_OrgBatchID','f_OrgName']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.orgbatch_id',
                   'orgbatch_id.organism_id.organism_name',
                   ]
 
#=================================================================================================
# WGS_CheckM - FastA CheckM
#=================================================================================================
class WGS_CheckM_Filter(Filterbase):
    f_OrgBatchID = CharFilter(field_name='orgbatch_id__orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
  
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def create_field_groups(self):
        self.group1 = [self.filters[name] for name in list(WGS_CheckM.HEADER_FIELDS.keys())]
        print(self.group1[0].label)

    class Meta:
        model = WGS_CheckM
        fields = ['f_OrgBatchID','f_OrgName']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.orgbatch_id',
                   'orgbatch_id.organism_id.organism_name',
                   ]

#=================================================================================================
# List of Genes
#=================================================================================================
## 
class Gene_Form(ModelForm):
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
 
class Gene_Filter(Filterbase):
   
    gene_type=ChoiceFilter(choices=[])
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["gene_type"].extra['choices']=[(obj.dict_value, obj) for obj in Dictionary.objects.filter(dict_class=Gene.Choice_Dictionary['gene_type'], astatus__gte=0)]
        
    class Meta:
        model=Gene
        fields=list(model.HEADER_FIELDS.keys())

#=================================================================================================
class AMRGenotype_Filter(Filterbase):
    f_OrgBatchID = CharFilter(field_name='orgbatch_id', lookup_expr='icontains',label="OrgBatch ID")
    f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
    #f_OrgName = CharFilter(field_name='orgbatch_id__organism_id__organism_name', lookup_expr='icontains',label="Organism")
    f_GeneCode = CharFilter(field_name='gene_id__gene_code', lookup_expr='icontains',label="Gene Code")
    f_GeneType = CharFilter(field_name='gene_id__gene_type', lookup_expr='icontains',label="Gene Type")
    f_GeneClass = CharFilter(field_name='gene_id__amr_class', lookup_expr='icontains',label="AMR Class")
    f_GeneSClass = CharFilter(field_name='gene_id__amr_subclass', lookup_expr='icontains',label="AMR SubClass")
    #gene_type=ChoiceFilter(choices=[])
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #self.filters["gene_type"].extra['choices']=[(obj.dict_value, obj) for obj in Dictionary.objects.filter(dict_class=Gene.Choice_Dictionary['gene_type'], astatus__gte=0)]
        
    class Meta:
        model=AMR_Genotype
        fields = ['f_OrgBatchID','f_OrgName','f_GeneCode','f_GeneType','f_GeneClass','f_GeneSClass']
        fields += list(model.HEADER_FIELDS.keys())
        exclude = ['orgbatch_id.orgbatch_id',
                    'orgbatch_id.organism_id.organism_name',
                    'gene_id.gene_code',
                    "gene_id.gene_type",
                    "gene_id.amr_class",
                    "gene_id.amr_subclass",
                   ]
        