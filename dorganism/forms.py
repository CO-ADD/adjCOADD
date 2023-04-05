from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField

from apputil.models import Dictionary, ApplicationUser
from .models import Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from adjcoadd.constants import *

# ----------------------------------------

class HiddenSimpleArrayField(forms.Field):
    widget = HiddenInput

    def clean(self, value):
        return value or []


#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):

    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    oxygen_pref=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['oxygen_pref'], astatus__gte=0), widget=forms.Select(attrs={'class': 'input-group'}), required=False,)
    risk_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['risk_group'], astatus__gte=0),widget=forms.Select(attrs={'class': 'input-group'}), required=False,)
    pathogen_group=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['pathogen_group'], astatus__gte=0),widget=forms.Select(attrs={'class': 'input-group'}),required=False,)
    mta_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'],  astatus__gte=0),widget=forms.Select(attrs={'class': 'input-group'}),required=False, )
    lab_restriction = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['lab_restriction']), widget=forms.Select(attrs={'class': 'input-group'}),required=False, )
    
    res_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    gen_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False,)
   
    def __init__(self, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        self.fields['strain_type'].widget = forms.SelectMultiple(choices= Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc=False),)
        self.fields['strain_type'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true',})
        self.fields['strain_panel'].widget = forms.SelectMultiple(choices= [])# Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['strain_panel'].widget.attrs.update({'class': 'form-select', 'size':'5', 'multiple': 'true'})
        self.fields['oxygen_pref'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['oxygen_pref'], astatus__gte=0)]
        self.fields['risk_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['risk_group'], astatus__gte=0)]
        self.fields['pathogen_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['pathogen_group'], astatus__gte=0)]
        self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'], astatus__gte=0)]
        self.fields['lab_restriction'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['lab_restriction'], astatus__gte=0)]
        self.create_field_groups()
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    def clean_organism_name(self):       
        data=self.cleaned_data['organism_name']
        data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
        return data

    def create_field_groups(self):
        self.group1 = [self[name] for name in ("strain_ids", "strain_code", "strain_notes", "strain_type", "strain_panel", "strain_origin", "strain_identification" )]
        self.group2 = [self[name] for name in ('res_property','gen_property','sequence_link','oxygen_pref','mta_status','mta_document', 'source', 'source_code')]
        self.group3 = [self[name] for name in ('risk_group','pathogen_group','lab_restriction','biologist','tax_id')]    
    
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
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['org_class'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class'], astatus__gte=0)]
        self.fields['division'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division'], astatus__gte=0)]

    class Meta:
        model =Taxonomy
        exclude=['urlname']
        fields=["organism_name","other_names", "code", "org_class", "tax_id", "parent_tax_id", "tax_rank", "division", "lineage" ]


#========================================Batch Form================================================================
class Batch_form(forms.ModelForm):
    # organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    qc_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0),required=False,)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    qc_record=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    
    def __init__(self, *args, **kwargs):
        super(Batch_form, self).__init__(*args, **kwargs)
        self.fields['qc_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0)]
        
      
              
    # def clean_organism_id(self):       
    #     data=self.cleaned_data['organism_id']
    #     try:
    #         data=get_object_or_404(Organism, organism_id=self.organism_id_str)#self.organism_name
    #     except Exception as err:
    #         print(err)
    #     return data
    
    

    class Meta:
        model =Organism_Batch
        exclude=['orgbatch_id', 'stock_level', 'organism_id']
        
# ---------------------------------------------------------------------------------------------
class Batchupdate_form(forms.ModelForm):
    qc_status = forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0),required=False,)
    orgbatch_id = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    qc_record=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    stock_level = HiddenSimpleArrayField(required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['qc_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Batch.Choice_Dictionary['qc_status'], astatus__gte=0)]
        self.create_field_groups()

    def create_field_groups(self):
        self.group_updatebatch = [self[name] for name in ("batch_id", "stock_date", "stock_level", "qc_status", "qc_record", "batch_notes", "biologist" )]
        
    class Meta:
        model =Organism_Batch
        fields=list(model.HEADER_FIELDS.keys())
        fields+=['orgbatch_id']
        exclude=['stock_level']


# ===============================Stock Form-------------------------------
class Stock_createform(forms.ModelForm):

    field_order = ['stock_type', 'n_created', 'n_left', 'stock_date', 'stock_note', 'passage_notes', 'location_freezer', 'location_rack', 'location_column', 'location_slot', 'biologist']

    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))
    # orgbatch_id=forms.ModelChoiceField(queryset=Organism_Batch.objects.filter(astatus__gte=0),widget=forms.HiddenInput())
    stock_type=forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=OrgBatch_Stock.Choice_Dictionary['stock_type'], astatus__gte=0), 
                                    widget=forms.Select(attrs={'class':'form-select', 'readonly':False}))
    passage_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    stock_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)

    
    def __init__(self, *args, **kwargs):
        # self.orgbatch_id = kwargs.pop('initial', None).get('orgbatch_id') if kwargs.get('initial') else None
        super().__init__(*args, **kwargs)
        self.fields['stock_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=OrgBatch_Stock.Choice_Dictionary['stock_type'], astatus__gte=0)]
        # self.fields['orgbatch_id'].choices=[(self.orgbatch_id, self.orgbatch_id)]

    class Meta:
        model =OrgBatch_Stock
        fields=list(model.HEADER_FIELDS.keys())
        exclude=['orgbatch_id', 'stock_type', 'stock_date','n_created']

class Stock_form(Stock_createform):
    # passage_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    class Meta:
        model =OrgBatch_Stock
        fields="__all__"

# ===============================Culture Form-------------------------------
class Culture_form(forms.ModelForm):
    # organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_type'], astatus__gte=0), 
                                    widget=forms.Select(attrs={'class':'form-select', 'readonly':False}), required=False,)
    culture_source= forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_source'], astatus__gte=0), 
                                    widget=forms.Select(attrs={'class':'form-select', 'readonly':False}), required=False,)

    # 
    def __init__(self, *args, **kwargs):

        super(Culture_form, self).__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_type'], astatus__gte=0)]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_source'], astatus__gte=0)]
     
              
    # def clean_organism_id(self):       
    #     data=self.cleaned_data['organism_id']
    #     data=get_object_or_404(Organism, organism_id=self.organism_id)#self.organism_name
    #     return data

    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys())
        exclude=['culture_type', 'culture_source'] 

# ---------------------------------------------------------------------------------------------
class Cultureupdate_form(forms.ModelForm):
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_type'], astatus__gte=0), 
                                    widget=forms.Select(attrs={'class':'form-select', 'width':'fit-content','disabled': 'disabled'}), required=False,)
    culture_source= forms.ModelChoiceField(queryset=Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_source'], astatus__gte=0), 
                                    widget=forms.Select(attrs={'class':'form-select',  'width':'fit-content','disabled': 'disabled'}), required=False,)
 
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_type'], astatus__gte=0)]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.objects.filter(dict_class=Organism_Culture.Choice_Dictionary['culture_source'], astatus__gte=0)]
        self.create_field_groups()

    def create_field_groups(self):
        self.group_updatebatch = [self[name] for name in ("culture_type", "culture_source", "media", "atmosphere", "temperature", "labware","culture_notes", "biologist" )]

        
    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys()) 
        exclude=['culture_type', 'culture_source',]

