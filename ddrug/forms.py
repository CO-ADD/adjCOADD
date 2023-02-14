from django import forms
from django.core.exceptions import ValidationError
from apputil.utils import get_DictonaryChoices_byDictClass
from django.shortcuts import get_object_or_404

from apputil.models import Dictionary, ApplicationUser
from .models import Drug
from adjcoadd.constants import *


#========================================Drug Form================================================================
class Drug_form(forms.ModelForm):
    drug_type = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Drug.Choice_Dictionary['drug_type']), widget=forms.Select(attrs={'class':'form-select'}))
    max_phase = forms.ChoiceField(choices= get_DictonaryChoices_byDictClass(Dictionary, Drug.Choice_Dictionary['max_phase'], ' | '), widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    class Meta:
        model =Drug
        fields='__all__'
        exclude=['ffp2']
