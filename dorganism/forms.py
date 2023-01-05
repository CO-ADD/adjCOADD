from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from apputil.utils import get_DictonaryChoices_byDictClass
from django.shortcuts import get_object_or_404

from apputil.models import Dictionary, ApplicationUser
from .models import Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from adjcoadd.constants import *

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
   
    def __init__(self, user, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        user=user
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.fields['strain_type'].widget = forms.SelectMultiple(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['strain_type'], ' | '))
        self.fields['strain_type'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['strain_panel'].widget = forms.SelectMultiple(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['strain_panel'], ' | '))
        self.fields['strain_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.initial['biologist']= get_object_or_404(ApplicationUser, name=user)
              
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
    org_class = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class']), widget=forms.Select(attrs={'class':'form-select'}))
    division = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    class Meta:
        model =Taxonomy
        exclude=['urlname']


#========================================Batch Form================================================================
class Batch_form(forms.ModelForm):
    # Batch_form_fields=ORGANISM_BATCH_FIELDs.keys()
    organism_id=forms.ModelChoiceField(queryset=Organism.objects.all(), widget=forms.HiddenInput(),required=False,)
   
    def __init__(self, user, organism_id=None, *args, **kwargs):
        self.organism_id=organism_id
        user=user
        super(Batch_form, self).__init__(*args, **kwargs)
        self.initial['biologist']=get_object_or_404(ApplicationUser, name=user)
              
    def clean_organism_id(self):       
        data=self.cleaned_data['organism_id']
        # organism=get_object_or_404(Taxonomy, organism_name=self.organism_name)
        data=get_object_or_404(Organism, organism_id=self.organism_id)#self.organism_name
        return data

    class Meta:
        model =Organism_Batch
        # fields=ORGANISM_BATCH_FIELDs.keys()
        # +["organism_id"]
        exclude=['orgbatch_id', 'stock_level']

class Batchupdate_form(forms.ModelForm):
    class Meta:
        model =Organism_Batch
        fields=ORGANISM_BATCH_FIELDs.keys()
        exclude=['stock_level']




# ===============================Stock Form-------------------------------
class Stock_form(forms.ModelForm):  
    class Meta:
        model =OrgBatch_Stock
        fields=ORGANISM_STOCK_FIELDs.keys()

# ===============================Culture Form-------------------------------
class Culture_form(forms.ModelForm):
    organism_id=forms.ModelChoiceField(queryset=Organism.objects.all(), widget=forms.HiddenInput(),required=False,)
    culture_type= forms.ChoiceField(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism_Culture.Choice_Dictionary['culture_type'], ' | '), widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    media_use= forms.ChoiceField(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism_Culture.Choice_Dictionary['media_use'], ' | '), widget=forms.Select(attrs={'class':'form-select'}), required=False,)

    def __init__(self, user, organism_id=None, *args, **kwargs):
        self.organism_id=organism_id
        user=user
        super(Culture_form, self).__init__(*args, **kwargs)
        self.initial['biologist']=get_object_or_404(ApplicationUser, name=user)
              
    def clean_organism_id(self):       
        data=self.cleaned_data['organism_id']
        # organism=get_object_or_404(Taxonomy, organism_name=self.organism_name)
        data=get_object_or_404(Organism, organism_id=self.organism_id)#self.organism_name
        return data

    class Meta:
        model =Organism_Culture
        fields=ORGANISM_CULTR_FIELDs

class Cultureupdate_form(forms.ModelForm):
    culture_type= forms.ChoiceField(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism_Culture.Choice_Dictionary['culture_type'], ' | '), widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    media_use= forms.ChoiceField(choices= get_DictonaryChoices_byDictClass(Dictionary, Organism_Culture.Choice_Dictionary['media_use'], ' | '), widget=forms.Select(attrs={'class':'form-select'}), required=False,)


    class Meta:
        model =Organism_Culture
        fields=ORGANISM_CULTR_FIELDs