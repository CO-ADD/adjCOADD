import django_filters
from rdkit import Chem
from django_rdkit.models import *

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.postgres.forms import SimpleArrayField
from django.shortcuts import get_object_or_404

from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import Filterbase
from .models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub
from adjcoadd.constants import *


#========================================Drug Form================================================================
class Drug_form(forms.ModelForm):
    drug_type = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['drug_type']), 
                                       widget=forms.Select(attrs={'class':'form-select'}), required=False)
    max_phase = forms.ChoiceField(choices= Dictionary.get_aschoices(Drug.Choice_Dictionary['max_phase'], showDesc=False), 
                                  widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    drug_codes= SimpleArrayField(forms.CharField(), required=False)
    drug_othernames = SimpleArrayField(forms.CharField(), required=False)
    drug_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    approval_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    drug_id=forms.CharField(widget=forms.HiddenInput(), required=False)
    # drug_class=forms.ChoiceField(choices=Dictionary.get_aschoices(Drug.Choice_Dictionary['drug_class'], showDesc=False), required=False)
    
    def __init__(self, *args, **kwargs):     
        super().__init__(*args, **kwargs)
        self.fields['drug_panel'].widget = forms.CheckboxSelectMultiple(choices= [])# Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['drug_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['drug_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['drug_type'], astatus__gte=0)]
        self.fields['max_phase'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['max_phase'], astatus__gte=0)]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    def create_field_groups(self):
        self.group1 = [self[name] for name in ("drug_othernames", "drug_codes", "drug_type", "drug_class", "drug_subclass", "drug_target", "drug_subtarget", "drug_panel","drug_note")]
        self.group2 = [self[name] for name in ('approval_note','admin_routes','application','n_compounds','chembl', 'drugbank', 'cas', 'pubchem', 'chemspider','unii', 'kegg', 'comptox', 'echa', 'chebi', 'uq_imb', 'vendor', 'vendor_catno')]
        self.group3 = [self[name] for name in ( 'moa', 'antimicro', 'antimicro_class','max_phase','mw','mf',)]


    class Meta:
        model =Drug
        fields='__all__'
        exclude=['ffp2', 'torsionbv', 'mfp2', 'smol']
    
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     for field in self.fields: 
    #         field.wiget.attrs['readonly'] = 'readonly'

    # def clean_smol(self):
    #     data=self.cleaned_data['smol']
    #     data=Chem.MolFromSmiles(data)
    #     print(data)
    #     return data

    # def clean_mfp2(self):
    #     data=self.cleaned_data['mfp2']
    #     data=Chem.AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(self.instance.smiles),radius=2, bitInfo={})
    #     print(data)
    #     return data

class Drug_updateform(Drug_form):
    # drug_name=forms.CharField(widget=forms.TextInput(attrs={'disabled': 'disabled'}),)

    def clean_unique_field(self):
        unique_field = self.cleaned_data['drug_name']
        # Exclude the current instance from the queryset to avoid the unique constraint conflict
        queryset =Drug.objects.exclude(pk=self.instance.pk)

        if queryset.filter(unique_field=unique_field).exists():
            raise forms.ValidationError("This unique_field value already exists.")
        
        return unique_field
    class Meta:
        model =Drug
        fields='__all__'
        exclude=['ffp2', 'torsionbv', 'mfp2', 'smol']


# -------------fitlerset Forms---------------------------------------------------------------

class Drug_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_name', lookup_expr='icontains')
    Drug_Type=django_filters.ChoiceFilter(field_name='drug_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Target=django_filters.CharFilter(field_name='drug_target', lookup_expr='icontains')
    Drug_Class=django_filters.CharFilter(field_name='drug_class', lookup_expr='icontains')
    Antimicro=django_filters.CharFilter(field_name='antimicro', lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Drug_Type"].extra['choices']=[(obj.dict_value, obj.__repr__()) for obj in Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['drug_type'], astatus__gte=0)]
        self.filters['Drug_Name'].label='Drug Name'
        self.filters['Drug_Type'].label='Drug Type'
        self.filters['Target'].label='Drug Target'
        self.filters['Drug_Class'].label='Drug Class'
        self.filters['Antimicro'].label='Antimicro'
    
    class Meta:
        model=Drug
        fields=['Drug_Name', 'Drug_Type', 'Target', 'Drug_Class', 'Antimicro']



class Vitekcard_filter(Filterbase):
    card_barcode = django_filters.CharFilter(lookup_expr='icontains')
    class Meta:
        model=VITEK_Card
        fields=['card_barcode']


class Vitekast_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_id__drug_name', lookup_expr='icontains')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['Drug_Name'].label='Drug Name'
    class Meta:
        model=VITEK_AST
        fields=['Drug_Name']


class MIC_COADDfilter(Filterbase):
    mic = django_filters.CharFilter(lookup_expr='icontains', label="MIC")
    class Meta:
        model=MIC_COADD
        fields=['mic']

class MIC_Pubfilter(Filterbase):
    mic = django_filters.CharFilter(lookup_expr='icontains', label="MIC")
    class Meta:
        model=MIC_Pub
        fields=['mic']