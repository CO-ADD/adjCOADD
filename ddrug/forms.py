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
    drug_type = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['drug_type']), widget=forms.Select(attrs={'class':'form-select'}))
    max_phase = forms.ChoiceField(choices= Dictionary.get_aschoices(Drug.Choice_Dictionary['max_phase'], showDesc=False), widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    drug_codes= SimpleArrayField(forms.CharField(), required=False)
    drug_othernames = SimpleArrayField(forms.CharField(), required=False)
    drug_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)

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

