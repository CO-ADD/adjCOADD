from django import forms
from django.core.paginator import Paginator
from django.forms import ModelForm
from dorganism.utils import querysetToChoiseList_Dictionary
from apputil.models import Dictionary, ApplicationUser
from .models import Organism, Taxonomy
from django.shortcuts import get_object_or_404
from django.core.exceptions import ValidationError
from django.contrib.postgres.forms import SimpleArrayField

#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):
 
    strain_ids= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    strain_code=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    strain_tissue=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    oxygen_pref=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['oxygen_pref']), widget=forms.Select(attrs={'class':'form-select'}))
    risk_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['risk_group']), widget=forms.Select(attrs={'class':'form-select'}))
    pathogen_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['pathogen_group']), widget=forms.Select(attrs={'class':'form-select'}))
    mta_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status']), widget=forms.Select(attrs={'class':'form-select'}))
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all())
    strain_panel=SimpleArrayField(forms.MultipleChoiceField(required=False, choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_type'])))
   
    def __init__(self, user, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        user=user
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.fields['strain_type'].widget = forms.SelectMultiple(choices= querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_type']))
        self.fields['strain_type'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['strain_panel'].widget = forms.SelectMultiple(choices= querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_panel']))
        self.fields['strain_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.initial['biologist']= ApplicationUser.objects.filter(username=user)[0]
              
    def clean_organism_name(self):       
        data=self.cleaned_data['organism_name']
        data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
        return data
            
    class Meta:
        model=Organism
        exclude = ['organism_id']

#=======================================Organism update Form=============================================================
class UpdateOrganism_form(CreateOrganism_form):       
    
    class Meta:
        model=Organism
        exclude = ['organism_id']

#========================================Taxonomy Form================================================================
class Taxonomy_form(forms.ModelForm):
    class Meta:
        model =Taxonomy
        fields='__all__'
