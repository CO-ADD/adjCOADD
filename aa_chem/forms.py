from .models import Organisms
from django import forms
from django.forms import ModelForm


class CreateNewOrgForm(ModelForm):
    Oxygen_Pref=forms.ChoiceField()
    Risk_Group=forms.ChoiceField()
    Pathogen=forms.ChoiceField()
   
   
    def __init__(self, choice1, choice2, choice3, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Oxygen_Pref_choices=choice1
        self.Risk_Group_choices=choice2
        self.Pathogen_choices=choice3
        
    #     # self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.testList)
        self.fields['Oxygen_Pref'].choices=self.Oxygen_Pref_choices
        self.fields['Risk_Group'].choices=self.Risk_Group_choices
        self.fields['Pathogen'].choices=self.Pathogen_choices

    class Meta:
        model=Organisms
        exclude = ['Strain_Type', 'Organism_Name']