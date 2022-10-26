from .models import Organisms, Taxonomy
from django import forms
from django.forms import ModelForm
from django.core.paginator import Paginator


class CreateNewOrgForm(ModelForm):
    Oxygen_Pref=forms.ChoiceField()
    Risk_Group=forms.ChoiceField()
    Pathogen=forms.ChoiceField()
    MTA_Status = forms.ChoiceField()
    Biol_Approval = forms.ChoiceField()
    # Organism_Name=forms.ModelChoiceField(queryset=Paginator(Taxonomy.objects.all(), 10) )
   
   
    def __init__(self, strain_type_choice, choice1, choice2, choice3, choice4, choice5, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.testList=strain_type_choice
        self.fields['Strain_Type'].widget = forms.CheckboxSelectMultiple(choices=self.testList)
        
        self.fields['Oxygen_Pref'].choices=choice1
        self.fields['Risk_Group'].choices=choice2
        self.fields['Pathogen'].choices=choice3
        self.fields['MTA_Status'].choices=choice4
        self.fields['Biol_Approval'].choices=choice5


    class Meta:
        model=Organisms
        exclude = ['Organism_Name', 'Organism_ID']

class UpdateNewOrgForm(ModelForm):
    
    class Meta:
        exclude = [ 'Organism_ID']
    
    # def get_context_data(self, **kwargs):
    #     context = super(UpdateNewOrgForm, self).get_context_data(**kwargs)
    #     comments = context['post'].comment_set.all()
    #     paginator = Paginator(comments, per_page=50)
    #     page_number = 1  # get it from query sting or whatever the way you want
    #     page = paginator.page(page_number)
    #     context['comments'] = page
        # return context
    


