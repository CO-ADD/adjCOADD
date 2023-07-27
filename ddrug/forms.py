import django_filters
from rdkit import Chem
from django_rdkit.models import *

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.postgres.forms import SimpleArrayField
from django.shortcuts import get_object_or_404

from django.contrib.postgres.search import TrigramSimilarity
from django.db.models.functions import Greatest

from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import Filterbase
from .models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
from adjcoadd.constants import *


#========================================Drug Form================================================================
class Drug_form(forms.ModelForm):
    drug_type = forms.ModelChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, queryset=Dictionary.objects.all())
    max_phase = forms.ChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, choices= [], )
    drug_codes= SimpleArrayField(forms.CharField(), required=False)
    drug_othernames = SimpleArrayField(forms.CharField(), required=False)
    drug_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    approval_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    drug_id=forms.CharField(widget=forms.HiddenInput(), required=False)
    smol=forms.CharField(widget=forms.TextInput(),)
   


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['drug_panel'].widget = forms.CheckboxSelectMultiple(choices= [])# Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['drug_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['drug_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['drug_type'])]
        self.fields['max_phase'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['max_phase'])]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    def create_field_groups(self):
        self.group1 = [self[name] for name in ("drug_othernames", "drug_codes", "drug_type", 'n_compounds',"drug_panel",'mw','mf',"drug_note", )]
        self.group2 = [self[name] for name in ("drug_class", "drug_subclass", "drug_target", "drug_subtarget", 'moa', 'antimicro', 'antimicro_class','max_phase','approval_note','admin_routes','application',)]
        self.group3 = [self[name] for name in ('chembl', 'drugbank', 'cas', 'pubchem', 'chemspider','unii', 'kegg', 'comptox', 'echa', 'chebi', 'uq_imb', 'vendor', 'vendor_catno')]

    class Meta:
        model =Drug
        fields='__all__'
        exclude=['ffp2', 'torsionbv', 'mfp2', ]
       
    

    def clean_smol(self):
        data=self.cleaned_data['smol']
        
        if data:
            data=Chem.MolFromSmiles(data)
        else:
            self.add_error('smol', 'Provide smol value, currently is None')
        
        return data

    # def clean_mfp2(self):
    #     data=self.cleaned_data['mfp2']
    #     data=Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(self.instance.smiles),radius=2, bitInfo={})
    #     print(data)
    #     return data


# -------------fitlerset Forms---------------------------------------------------------------

class Drug_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_name', lookup_expr='icontains')
    Drug_Type=django_filters.ChoiceFilter(field_name='drug_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Target=django_filters.CharFilter(field_name='drug_target', lookup_expr='icontains')
    Drug_Class=django_filters.CharFilter(field_name='drug_class', lookup_expr='icontains')
    Antimicro=django_filters.CharFilter(field_name='antimicro', lookup_expr='icontains')
    Other_Name = django_filters.CharFilter(field_name='drug_othernames', lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Drug_Type"].extra['choices']=[(obj.dict_value, obj.__repr__()) for obj in Dictionary.get_filterobj(Drug.Choice_Dictionary['drug_type'])]
        self.filters['Drug_Name'].label='Drug Name'
        self.filters['Drug_Type'].label='Drug Type'
        self.filters['Target'].label='Drug Target'
        self.filters['Drug_Class'].label='Drug Class'
        self.filters['Antimicro'].label='Antimicro'
    
    class Meta:
        model=Drug
        fields=['Drug_Name', 'Drug_Type', 'Target', 'Drug_Class', 'Antimicro', 'Other_Name']
        # fields=list(model.HEADER_FIELDS.keys())



class Vitekcard_filter(Filterbase):
    card_barcode = django_filters.CharFilter(lookup_expr='icontains')
    class Meta:
        model=VITEK_Card
        fields=['card_barcode']


class Vitekast_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_id__drug_name', lookup_expr='icontains')
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.filters['Drug_Name'].label='Drug Name'
    class Meta:
        model=VITEK_AST
        fields=['Drug_Name']

class VitekID_filter(Filterbase):
    
    class Meta:
        model=VITEK_ID
        fields=list(model.HEADER_FIELDS.keys())

class MIC_COADDfilter(Filterbase):
    mic = django_filters.CharFilter(lookup_expr='icontains', label="MIC")
    # orgbatch_id__organism_id__organism_name = django_filters.CharFilter(lookup_expr='icontains')

    def filter_all_fields(self, queryset, name, value):
        if value:
            # print(f"filtr fields is {self._meta.model._meta.fields})")
            fields=[f.name for f in self._meta.model._meta.fields]
            value=value         
            similarity = Greatest(
                TrigramSimilarity('drug_id__drug_name', value),
                TrigramSimilarity('mic', value),
                TrigramSimilarity('orgbatch_id__organism_id__organism_name', value),
                TrigramSimilarity('mic_unit', value),
                TrigramSimilarity('mic_type__dict_value', value),
                TrigramSimilarity('bp_profile', value),
                TrigramSimilarity('bp_source', value),
                TrigramSimilarity('run_id', value),
                TrigramSimilarity('testplate_id', value),
                TrigramSimilarity('testwell_id', value),
                TrigramSimilarity('plate_size__dict_value', value),
                TrigramSimilarity('plate_material__dict_value', value),
                TrigramSimilarity('media', value),
                TrigramSimilarity('dye', value),


               )
            queryset=queryset.annotate(similarity=similarity)#(similarity=TrigramSimilarity('drug_id__drug_name', value),)

            return queryset.filter(similarity__gt=0.1)#(q_object)
        return queryset
    
    def filter_all_fields_deep(self, queryset, name, value):
        queryset=self.filter_all_fields(queryset, name, value)
        return queryset
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['drug_id__drug_name'].label='Drug Name'
        self.filters['orgbatch_id__organism_id__gen_property'].label='Gen Property'
        self.filters['mic'].label='MIC'

    class Meta:
        model=MIC_COADD
        fields= ["drug_id__drug_name", "orgbatch_id__organism_id__gen_property","mic"]

class MIC_Pubfilter(Filterbase):
    mic = django_filters.CharFilter(lookup_expr='icontains', label="MIC")
    def filter_all_fields(self, queryset, name, value):
        if value:
            # print(f"filtr fields is {self._meta.model._meta.fields})")
            fields=[f.name for f in self._meta.model._meta.fields]
            value=value         
            similarity = Greatest(
                TrigramSimilarity('organism_id__organism_name', value),
                TrigramSimilarity('drug_id__drug_name', value),
                TrigramSimilarity('mic', value),
                TrigramSimilarity('mic_unit', value),
                TrigramSimilarity('zone_diameter', value),
                TrigramSimilarity('mic_type__dict_value', value),
                TrigramSimilarity('source', value),
                TrigramSimilarity('bp_profile', value),
                TrigramSimilarity('bp_source', value),

               )
            queryset=queryset.annotate(similarity=similarity)#(similarity=TrigramSimilarity('drug_id__drug_name', value),)

            return queryset.filter(similarity__gt=0.1)#(q_object)
        return queryset
    
    def filter_all_fields_deep(self, queryset, name, value):
        queryset=self.filter_all_fields(queryset, name, value)
        return queryset

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['drug_id__drug_name'].label='Drug Name'
        self.filters['organism_id__gen_property'].label='Gen Property'
        self.filters['mic'].label='MIC'
    
    class Meta:
        model=MIC_Pub
        fields=["drug_id__drug_name", "organism_id__gen_property", "mic"]


class Breakpointfilter(Filterbase):
    drug_name = django_filters.CharFilter(field_name='drug_id.drug_name', lookup_expr='icontains', label="Drug")
    class Meta:
        model=Breakpoint
        fields = ['drug_name',]
        fields += list(model.HEADER_FIELDS.keys()) 
        exclude=['drug_id.drug_name']
 
