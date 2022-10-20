from .models import Organisms
from django import forms


class CreateNewOrgForm(forms.ModelForm):
    
    Strain_Type=forms.MultipleChoiceField(choices=Organisms.Choice_Dictionaries)
   
    def __init__(self, list_test, *args, **kwargs):
        self.testList=list_test
        super(CreateNewOrgForm, self).__init__(*args, **kwargs)
        # self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple()
        self.fields['Strain_Type'].choices=self.testList

    class Meta:
        model=Organisms
        fields='__all__'