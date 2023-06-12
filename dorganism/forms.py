import django_filters
from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField

from apputil.models import Dictionary, ApplicationUser, Image, Document
from apputil.utils.filters_base import Filterbase
from .models import Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from adjcoadd.constants import *

# ----------------------------------------

class HiddenSimpleArrayField(forms.Field):
    widget = HiddenInput

    def clean(self, value):
        return value or []


#=======================================Organism Create Form=============================================================
class CreateOrganism_form(ModelForm):

    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    oxygen_pref=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    risk_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    pathogen_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False,queryset=Dictionary.objects.all())
    mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    lab_restriction = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    
    res_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    gen_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=False,)
    
   
    def __init__(self, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        super(CreateOrganism_form, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['strain_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc=False),)
        self.fields['strain_type'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true',})
        self.fields['strain_panel'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc=False),)
        self.fields['strain_panel'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true'})
        self.fields['oxygen_pref'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['oxygen_pref'])]
        self.fields['risk_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['risk_group'])]
        self.fields['pathogen_group'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['pathogen_group'])]
        self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['mta_status'])]
        self.fields['lab_restriction'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism.Choice_Dictionary['lab_restriction'])]
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
        if data:
            return data
        else:
            self.add_error('organism_name', "Found No Organism")
            raise ValidationError

    def create_field_groups(self):
        self.group1 = [self[name] for name in Organism.FORM_GROUPS['Group1']]
        self.group2 = [self[name] for name in Organism.FORM_GROUPS['Group2']]
        self.group3 = [self[name] for name in Organism.FORM_GROUPS['Group3']]
        self.group4 = [self[name] for name in Organism.FORM_GROUPS['Group4']] 
    
    class Meta:
        model=Organism
        exclude=['organism_id', 'assoc_images', 'assoc_documents'] 

#=======================================Organism update Form=============================================================
class UpdateOrganism_form(CreateOrganism_form):   
    
    class Meta:
        model=Organism
        exclude=['organism_id', 'assoc_images', 'assoc_documents'] 
   
#========================================Taxonomy Form================================================================
class Taxonomy_form(forms.ModelForm):
    org_class = forms.ModelChoiceField(widget=forms.Select(attrs={'class':''}), queryset=Dictionary.objects.all(), )
    division = forms.ModelChoiceField(queryset=Dictionary.objects.all())
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['org_class'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Taxonomy.Choice_Dictionary['org_class'])]
        self.fields['division'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Taxonomy.Choice_Dictionary['division'])]

    class Meta:
        model =Taxonomy
        exclude=['urlname']
        fields=["organism_name","other_names", "code", "org_class", "tax_id", "parent_tax_id", "tax_rank", "division", "lineage" ]


#========================================Batch Form================================================================
class Batch_form(forms.ModelForm):
    # organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    qc_record=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    batch_id=forms.CharField(widget=forms.TextInput(attrs={'maxlength': '5', 'default':'optional input'}), help_text='**Optional with up to 5 characters', required=False)

    def __init__(self, *args, **kwargs):
        super(Batch_form, self).__init__(*args, **kwargs)
        self.fields['qc_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['qc_status'])]
    

    class Meta:
        model =Organism_Batch
        exclude=['orgbatch_id', 'stock_level', 'organism_id']
        
#========================================Batch Update Form================================================================
class Batchupdate_form(forms.ModelForm):
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    orgbatch_id = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    qc_record=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    stock_level = SimpleArrayField(forms.IntegerField(), delimiter=';', disabled=True)

    def __init__(self, *args, **kwargs):
        
        super().__init__(*args, **kwargs)
        instance=kwargs.get('instance')
        if instance and instance.stock_level:
            self.fields['stock_level'].initial=instance.stock_level
        self.fields['qc_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['qc_status'])]
        self.create_field_groups()

    def create_field_groups(self):
        self.group_updatebatch = [self[name] for name in ("batch_id", "stock_date", "stock_level", "qc_status", "qc_record", "batch_notes", "biologist" )]
        
    class Meta:
        model =Organism_Batch
        fields=list(model.HEADER_FIELDS.keys())
        fields+=['orgbatch_id']
        exclude=['stock_level']


# ===============================Stock Create Form-------------------------------
class Stock_createform(forms.ModelForm):

    field_order = ['stock_type', 'n_created', 'n_left', 'stock_date', 'stock_note', 'passage_notes', 'location_freezer', 'location_rack', 'location_column', 'location_slot', 'biologist']

    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))
    # orgbatch_id=forms.ModelChoiceField(queryset=Organism_Batch.objects.filter(astatus__gte=0),widget=forms.HiddenInput())
    stock_type=forms.ModelChoiceField(widget=forms.Select(attrs={'class':'', 'readonly':False}),queryset=Dictionary.objects.all(),)
    passage_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    stock_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)

    
    def __init__(self, *args, **kwargs):
        # self.orgbatch_id = kwargs.pop('initial', None).get('orgbatch_id') if kwargs.get('initial') else None
        super().__init__(*args, **kwargs)
        self.fields['stock_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(OrgBatch_Stock.Choice_Dictionary['stock_type'])]
        # self.fields['orgbatch_id'].choices=[(self.orgbatch_id, self.orgbatch_id)]

    class Meta:
        model =OrgBatch_Stock
        fields=list(model.HEADER_FIELDS.keys())
        exclude=['orgbatch_id', 'stock_type', 'stock_date','n_created']

#======================================== Stock Form================================================================
class Stock_form(Stock_createform):
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['stock_date'].initial=self.instance.stock_date
        self.fields['n_created'].initial=self.instance.n_created

    class Meta:
        model =OrgBatch_Stock
        fields="__all__"

# =============================== Culture Form-------------------------------
class Culture_form(forms.ModelForm):
    # organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(widget=forms.Select(attrs={'readonly':False}), required=False,queryset=Dictionary.objects.all(),)
    culture_source= forms.ModelChoiceField(widget=forms.Select(attrs={'readonly':False}), required=False,queryset=Dictionary.objects.all(),)

    # 
    def __init__(self, *args, **kwargs):

        super(Culture_form, self).__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_type'])]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_source'])]
     
              
    # def clean_organism_id(self):       
    #     data=self.cleaned_data['organism_id']
    #     data=get_object_or_404(Organism, organism_id=self.organism_id)#self.organism_name
    #     return data

    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys())
        # exclude=['culture_type', 'culture_source'] 

# =============================== Culture Update Form-------------------------------
class Cultureupdate_form(forms.ModelForm):
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(widget=forms.Select(attrs={'class':'', 'width':'fit-content','disabled': 'disabled'}), required=False,queryset=Dictionary.objects.all(),)
    culture_source= forms.ModelChoiceField(widget=forms.Select(attrs={'class':'',  'width':'fit-content','disabled': 'disabled'}), required=False,queryset=Dictionary.objects.all(),)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_type'])]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_source'])]
        self.create_field_groups()

    def create_field_groups(self):
        self.group_updatebatch = [self[name] for name in ("culture_type", "culture_source", "media", "atmosphere", "temperature", "labware","culture_notes", "biologist" )]

        
    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys()) 
        exclude=['culture_type', 'culture_source',]

# --Filterset Forms--
## Taxonomy
class Taxonomyfilter(Filterbase):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    lineage = django_filters.CharFilter(lookup_expr='icontains')
    org_class=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class']))
    division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['organism_name'].label='Organism Name'
        self.filters['code'].label='Code'
        self.filters['lineage'].label='Lineage'
        self.filters['tax_rank'].label='Rank'
        self.filters['division'].label='Division'
        self.filters['org_class'].label='Class'
        self.filters['tax_id'].label='Tax ID'

    class Meta:
        model=Taxonomy
        fields=['organism_name', 'code', 'lineage', 'tax_rank','division', 'org_class', 'tax_id', ]

## Organism
class Organismfilter(Filterbase):
    
    ID=django_filters.CharFilter(field_name='organism_id', lookup_expr='icontains')
    Name = django_filters.CharFilter(field_name='organism_name__organism_name', lookup_expr='icontains')
    Class=django_filters.ChoiceFilter(field_name='organism_name__org_class__dict_value',  
                                      widget=forms.RadioSelect, 
                                      choices=(("GN","GN"),("GP","GP"),("FG","FG"),("MB","MB")))
    Strain=django_filters.CharFilter(field_name='strain_ids', lookup_expr='icontains')
    Notes=django_filters.CharFilter(field_name='strain_notes', lookup_expr='icontains')
    Type=django_filters.MultipleChoiceFilter(field_name='strain_type', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    MTA=django_filters.ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'], astatus__gte=0))
    Panel=django_filters.MultipleChoiceFilter(field_name='strain_panel',method='multichoices_filter', choices=Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc = False))
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Type"].extra["choices"]=Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc = False)
        for i in self.filters:
            self.filters[i].label=i
    
    
    class Meta:
        model=Organism
        fields=[ 'Class', 'ID', 'Name','Strain',  'Notes', 'Type', 'MTA', 'Panel', ]
       


## Batch
class Batchfilter(Filterbase):
    Stock_Date=django_filters.IsoDateTimeFilter(field_name='stock_date')
    class Meta:
        model=Organism_Batch
        fields= ["stock_date",  "biologist"]


