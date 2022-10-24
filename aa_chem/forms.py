from .models import Organisms
from django import forms
from django.forms import ModelForm



class CreateNewOrgForm(ModelForm):
    
    # Strain_Type=forms.MultipleChoiceField()
   
    # def __init__(self, list_test, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.testList=list_test
        
    #     # self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.testList)
    #     self.fields['Strain_Type'].choices=self.testList

    class Meta:
        model=Organisms
        exclude = ['Strain_Type', 'Organism_Name']