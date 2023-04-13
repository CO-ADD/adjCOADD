from rdkit import Chem
from django_rdkit.models import *

from django import forms
from django.core.exceptions import ValidationError
from django.contrib.postgres.forms import SimpleArrayField
from django.shortcuts import get_object_or_404

from apputil.models import Dictionary, ApplicationUser
from .models import Drug
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
        self.group1 = [self[name] for name in ("drug_othernames", "drug_codes", "drug_type", "drug_class", "drug_subclass", "drug_target", "drug_subtarget", "drug_panel",)]
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

