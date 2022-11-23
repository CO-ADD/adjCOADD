from django import forms
from django.core.paginator import Paginator
from django.forms import ModelForm
from dorganism.utils import querysetToChoiseList_Dictionary
from apputil.models import Dictionary, ApplicationUser

from .models import Organism, Taxonomy
from django.shortcuts import get_object_or_404
from django.core.exceptions import ValidationError

#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):

    # Organism_Desc=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    strain_id= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    strain_code=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    # Strain_Desc= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    strain_tissue=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    oxygen_pref=forms.ChoiceField(choices=(('',''),('','')), widget=forms.Select(attrs={'class':'form-select special'},))
    risk_group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    pathogen_group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    mta_status = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    bio_approval = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all())
   
    def __init__(self, user, organism_name=None, *args, **kwargs): #Strain_Type_choices,
        self.organism_name=organism_name
        user=user
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        strain_type_choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_type']) # 
        # user=kwargs.pop("user")
        

        # print(type(user))
        self.initial['biologist']= ApplicationUser.objects.filter(username=user)[0]
       
        self.strainTypeChoices= strain_type_choices
        
        self.fields['strain_type'].widget = forms.SelectMultiple(choices=self.strainTypeChoices)
        self.fields['strain_type'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['oxygen_pref'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['oxygen_pref'])
        self.fields['risk_group'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['risk_group'])
        self.fields['pathogen_group'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['pathogen_group'])
        self.fields['mta_status'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['mta_status'])
        self.fields['bio_approval'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['bio_approval'])
           
    def clean_Organism_Name(self):       
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




    
    
    


