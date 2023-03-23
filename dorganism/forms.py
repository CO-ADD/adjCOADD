from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404

from apputil.models import Dictionary, ApplicationUser
from .models import Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from adjcoadd.constants import *

#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):
 
    strain_ids= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    strain_code=forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}), required=False,)
    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    oxygen_pref=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['oxygen_pref'], astatus__gte=0), widget=forms.Select(attrs={}), required=False,)
    risk_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['risk_group'], astatus__gte=0),widget=forms.Select(attrs={}), required=False,)
    pathogen_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['pathogen_group'], astatus__gte=0),widget=forms.Select(attrs={}),required=False,)
    mta_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'],  astatus__gte=0),widget=forms.Select(attrs={}),required=False, )
    lab_restriction = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['lab_restriction']), widget=forms.Select(attrs={}),required=False, )

    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False,)
   
    def __init__(self, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.fields['strain_type'].widget = forms.SelectMultiple(choices= Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc=False),)
        self.fields['strain_type'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true',})
        self.fields['strain_panel'].widget = forms.SelectMultiple(choices= [])# Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['strain_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['mta_status'].queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'], astatus__gte=0)
              
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
    org_class = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class'], astatus__gte=0), widget=forms.Select(attrs={'class':'form-select'}))
    division = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division'], astatus__gte=0))
    class Meta:
        model =Taxonomy
        exclude=['urlname']
        fields=["organism_name","other_names", "code", "org_class", "tax_id", "parent_tax_id", "tax_rank", "division", "lineage" ]


#========================================Batch Form================================================================
class Batch_form(forms.ModelForm):
    organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    qc_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0),required=False,)
    def __init__(self, organism_id_str=None, *args, **kwargs):
        self.organism_id_str=organism_id_str
        super(Batch_form, self).__init__(*args, **kwargs)
      
              
    def clean_organism_id(self):       
        data=self.cleaned_data['organism_id']
        try:
            data=get_object_or_404(Organism, organism_id=self.organism_id_str)#self.organism_name
        except Exception as err:
            print(err)
        return data

    class Meta:
        model =Organism_Batch
        exclude=['orgbatch_id', 'stock_level']
        
# ---------------------------------------------------------------------------------------------
class Batchupdate_form(forms.ModelForm):
    qc_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0),required=False,)
    orgbatch_id = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),)
    class Meta:
        model =Organism_Batch
        fields=model.HEADER_FIELDS.keys()
        exclude=['stock_level']


# ===============================Stock Form-------------------------------
class Stock_createform(forms.ModelForm):
    n_left_extra=forms.IntegerField(required=True)  
    class Meta:
        model =OrgBatch_Stock
        fields="__all__"

class Stock_form(forms.ModelForm):
    class Meta:
        model =OrgBatch_Stock
        fields="__all__"

# ===============================Culture Form-------------------------------
class Culture_form(forms.ModelForm):
    organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    # culture_type= forms.ChoiceField(choices= Dictionary.get_aschoices(Organism_Culture.Choice_Dictionary['culture_type'], showDesc=True), 
    #                                 widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    # media_use= forms.ChoiceField(choices= Dictionary.get_aschoices(Organism_Culture.Choice_Dictionary['media_use'], showDesc=True),
    # 
    def __init__(self, organism_id=None, *args, **kwargs):
        self.organism_id=organism_id
        super(Culture_form, self).__init__(*args, **kwargs)
        
              
    def clean_organism_id(self):       
        data=self.cleaned_data['organism_id']
        data=get_object_or_404(Organism, organism_id=self.organism_id)#self.organism_name
        return data

    class Meta:
        model =Organism_Culture
        fields=model.HEADER_FIELDS 

# ---------------------------------------------------------------------------------------------
class Cultureupdate_form(forms.ModelForm):
    # culture_type= forms.ChoiceField(choices= Dictionary.get_aschoices(Organism_Culture.Choice_Dictionary['culture_type'], showDesc=True), 
    #                                 widget=forms.Select(attrs={'class':'form-select'}), required=False,)
    # media_use= forms.ChoiceField(choices= Dictionary.get_aschoices(Organism_Culture.Choice_Dictionary['media_use'], showDesc=True),
    # 
    class Meta:
        model =Organism_Culture
        fields=model.HEADER_FIELDS 