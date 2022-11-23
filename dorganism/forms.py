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
    Strain_ID= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    Strain_Code=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    # Strain_Desc= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    Strain_Tissue=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    Oxygen_Pref=forms.ChoiceField(choices=(('',''),('','')), widget=forms.Select(attrs={'class':'form-select special'},))
    Risk_Group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Pathogen_Group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    MTA_Status = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Bio_Approval = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Organism_Name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    Biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all())
   
    def __init__(self, user, Organism_Name=None, *args, **kwargs): #Strain_Type_choices,
        self.Organism_Name=Organism_Name
        user=user
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        Strain_Type_choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['Strain_Type']) # 
        # user=kwargs.pop("user")
        

        # print(type(user))
        self.initial['Biologist']= ApplicationUser.objects.filter(username=user)[0]
       
        self.strainTypeChoices= Strain_Type_choices
        
        self.fields['Strain_Type'].widget = forms.SelectMultiple(choices=self.strainTypeChoices)
        self.fields['Strain_Type'].widget.attrs.update({'class': 'form-select', 'size':'3', 'multiple': 'true'})
        self.fields['Oxygen_Pref'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['Oxygen_Pref'])
        self.fields['Risk_Group'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['Risk_Group'])
        self.fields['Pathogen_Group'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['Pathogen_Group'])
        self.fields['MTA_Status'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['MTA_Status'])
        self.fields['Bio_Approval'].choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['Bio_Approval'])
           
    def clean_Organism_Name(self):       
        data=self.cleaned_data['Organism_Name']
        data=get_object_or_404(Taxonomy, Organism_Name=self.Organism_Name)
        return data
            
    class Meta:
        model=Organism
        exclude = ['Organism_ID']

#=======================================Organism update Form=============================================================
class UpdateOrganism_form(CreateOrganism_form):       
    
    class Meta:
        model=Organism
        exclude = ['Organism_ID']

#========================================Taxonomy Form================================================================
class Taxonomy_form(forms.ModelForm):
    class Meta:
        model =Taxonomy
        fields='__all__'




    
    
    


