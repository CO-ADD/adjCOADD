from django import forms
from django.core.paginator import Paginator
from django.forms import ModelForm
from aa_chem.utils import querysetToChoiseList_Dictionaries
from app.models import Dictionaries
from .models import Organisms, Taxonomy
from django.shortcuts import get_object_or_404
from django.core.exceptions import ValidationError

#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):
   
    Organism_Desc=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_ID= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Code=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Desc= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Notes= forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Strain_Tissue=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control'}), required=False,)
    Oxygen_Pref=forms.ChoiceField(choices=(('',''),('','')), widget=forms.Select(attrs={'class':'form-select', 'id':'Oxygen_Pref_choice'},))
    Risk_Group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Pathogen_Group=forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    MTA_Status = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Bio_Approval = forms.ChoiceField(choices=(('',''),('','')),  widget=forms.Select(attrs={'class':'form-select'}))
    Organism_Name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
   
   
    def __init__(self, Strain_Type_choices,  *args, **kwargs):
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.strainTypeChoices= Strain_Type_choices
        self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.strainTypeChoices)
        self.fields['Oxygen_Pref'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Oxygen_Pref'])
        self.fields['Risk_Group'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Risk_Group'])
        self.fields['Pathogen_Group'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Pathogen_Group'])
        self.fields['MTA_Status'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['MTA_Status'])
        self.fields['Bio_Approval'].choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Bio_Approval'])
        
    
    def get_object(self, Organism_Name):
        Organism =Taxonomy.objects.filter(Organism_Name=Organism_Name) #get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
        # print(Organism.Class.Dict_Value)
        
        self.fields['Organism_Name'].queryset=Organism
        print(f'Thisis from def get_object: {self.fields["Organism_Name"]}')
        return self.fields['Organism_Name']
       

            
    class Meta:
        model=Organisms
        exclude = ['Organism_ID']

#=======================================Organism update Form=============================================================
class UpdateOrganism_form(CreateOrganism_form):

    def clean_organismName(self, Organism_Name, original_class):
        print("start clean...")
        newOrganism=get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
        try:
            newOrganism_class=newOrganism.Class.Dict_Value
        except Exception as err:
            print(err)
            newOrganism_class=" "   
        print(newOrganism_class, original_class)
        if original_class != newOrganism_class:
            raise ValidationError("New Organism has different Class, which is not allowed!")
        
        self.fields['Organism_Name']=newOrganism
        # print(self.fields['Organism_Name'])
        # return self.fields['Organism_Name']
        
    
    class Meta:
        model=Organisms
        exclude = ['Organism_ID']

#========================================Taxonomy Form================================================================
class Taxonomy_form(forms.ModelForm):
    class Meta:
        model =Taxonomy
        fields='__all__'




    
    
    


