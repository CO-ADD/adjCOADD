import django_filters
from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField

# from apputil.models import Dictionary, ApplicationUser, Document
# from apputil.utils.filters_base import Filterbase
# from .models import Newmodel ...
# from adjcoadd.constants import *

# ----------------------------------------

class HiddenSimpleArrayField(forms.Field):
    widget = HiddenInput

    def clean(self, value):
        return value or []
class DateField(forms.Field):
    widget = forms.DateInput(attrs={'type': 'date'})
    def clean(self, value):
        return value

class NewModelimg_form(forms.ModelForm):
    pass
    # field_order = ['image_file','orgbatch_id','image_desc','image_source']
    # image_file = forms.ImageField(label='Select an image file',required=True)
    # image_type = forms.CharField(widget=forms.TextInput(attrs={'readonly':'readonly'}))
    # image_name = forms.CharField(widget=forms.TextInput(attrs={'readonly':'readonly'}))
 
    # def __init__(self, *args, org=None, **kwargs):
        
    #     super().__init__(*args, **kwargs)
    #     if org:
    #         pk = org
    #         organism = get_object_or_404(Organism, pk=pk)
    #         self.fields['orgbatch_id'].queryset = Organism_Batch.objects.filter(organism_id = organism.pk)
        
    # class Meta:
    #     model =OrgBatch_Image
    #     fields="__all__"
#=======================================New Model Create Form=============================================================
class CreateNewModel_form(forms.ModelForm):
    pass

    # strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    # # prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    # oxygen_pref=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    # risk_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    # pathogen_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False,queryset=Dictionary.objects.all())
    # mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    # lab_restriction = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    # res_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    # gen_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    # organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    # biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    # collect_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
   
    # def __init__(self, organism_name=None, *args, **kwargs): 
    #     self.organism_name=organism_name
    #     super(CreateOrganism_form, self).__init__(*args, **kwargs)
    #     for field_name in self.fields:
    #         self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
    #     self.fields['strain_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc=False),)
    #     self.fields['strain_type'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true',})
    #     self.fields['strain_panel'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
    #     self.fields['strain_panel'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true'})
    #     self.fields['oxygen_pref'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['oxygen_pref'])]
    #     self.fields['risk_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['risk_group'])]
    #     self.fields['pathogen_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['pathogen_group'])]
    #     self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['mta_status'])]
    #     self.fields['lab_restriction'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['lab_restriction'])]
    #     self.create_field_groups()

    #     for field in self.fields.values():
    #         if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
    #             # Add the 'group-input' class to the widget attrs
    #             attrs = field.widget.attrs
    #             attrs['class'] = attrs.get('class', '') + 'input-group'
    #             field.widget.attrs = attrs
    
    # def clean_organism_name(self):       
       
    #     data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
    
    #     if data:
    #         if str(data.org_class) in ORGANISM_CLASSES:
    #             return data
    #         else:
    #             self.add_error('organism_name', 'Create failed with invalid organism class')
                
    #     else:
    #         self.add_error('organism_name', "Found No Organism")

    #     return data            


    # def create_field_groups(self):
    #     self.group1 = [self[name] for name in Organism.FORM_GROUPS['Group1']]
    #     self.group2 = [self[name] for name in Organism.FORM_GROUPS['Group2']]
    #     self.group3 = [self[name] for name in Organism.FORM_GROUPS['Group3']]
    #     self.group4 = [self[name] for name in Organism.FORM_GROUPS['Group4']] 
    
    # class Meta:
    #     model=Organism
    #     exclude=['organism_id',  'assoc_documents'] 

#=======================================Organism update Form=============================================================
class UpdateNewModel_form(CreateNewModel_form): 
    pass
    
    # class Meta:
    #     model=Organism
    #     exclude=['organism_id', 'assoc_documents'] 
   


# --Filterset Forms--
## NewModelfilter
class NewModelfilter(Filterbase):
    pass
    # organism_name = django_filters.CharFilter(lookup_expr='icontains')
    # lineage = django_filters.CharFilter(lookup_expr='icontains')
    # org_class=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class']))
    # division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    
    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.filters['organism_name'].label='Organism Name'
    #     self.filters['code'].label='Code'
    #     self.filters['lineage'].label='Lineage'
    #     self.filters['tax_rank'].label='Rank'
    #     self.filters['division'].label='Division'
    #     self.filters['org_class'].label='Class'
    #     self.filters['tax_id'].label='Tax ID'

    # class Meta:
    #     model=Taxonomy
    #     fields=['organism_name', 'code', 'lineage', 'tax_rank','division', 'org_class', 'tax_id', ]

