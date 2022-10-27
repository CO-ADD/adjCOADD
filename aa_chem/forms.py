from .models import Organisms, Taxonomy
from app.models import Dictionaries
from django import forms
from django.forms import ModelForm
from django.core.paginator import Paginator
from aa_chem.utils import  querysetToChoiseList_Dictionaries

#=======================================Organism Create Form=============================================================
class CreateNewOrgForm(ModelForm):
    Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type'])
    Oxygen_Pref_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Oxygen_Pref']) 
    Risk_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Risk_Group'])
    MTA_Status_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['MTA_Status'])
    Biological_Approval_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Bio_Approval'])
    Pathogen_Group_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Pathogen_Group'])
   
    Organism_Desc=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_ID= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Code=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Desc= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Notes= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Tissue=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Oxygen_Pref=forms.ChoiceField(choices=Oxygen_Pref_choices, widget=forms.Select(attrs={'class':'form-select'}))
    Risk_Group=forms.ChoiceField(choices=Risk_Group_choices,  widget=forms.Select(attrs={'class':'form-select'}))
    Pathogen=forms.ChoiceField(choices=Pathogen_Group_choices,  widget=forms.Select(attrs={'class':'form-select'}))
    MTA_Status = forms.ChoiceField(choices=MTA_Status_choices,  widget=forms.Select(attrs={'class':'form-select'}))
    Biol_Approval = forms.ChoiceField(choices=Biological_Approval_choices,  widget=forms.Select(attrs={'class':'form-select'}))
    # Organism_Name=forms.ModelChoiceField(queryset=Paginator(Taxonomy.objects.all(), 10) )
   
   
    def __init__(self, Strain_Type_choices, *args, **kwargs):
        super(CreateNewOrgForm, self).__init__(*args, **kwargs)
        self.testList=Strain_Type_choices
        self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.testList)
        
            
    class Meta:
        Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type'])
        model=Organisms
        exclude = ['Organism_Name', 'Organism_ID']

       
#=======================================Organism update Form=============================================================
class UpdateNewOrgForm(CreateNewOrgForm):
    
    class Meta:
        model=Organisms
        exclude = ['Organism_Name', 'Organism_ID']





#========================================Taxonomy Form================================================================
class TaxonomyCreateForm(forms.ModelForm):
    class Meta:
        model =Taxonomy
        fields='__all__'
    
    
    


