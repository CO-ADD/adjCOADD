from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter

from apputil.models import Dictionary, ApplicationUser, Document
from apputil.utils.filters_base import Filterbase
from adjcoadd.constants import ORGANISM_CLASSES
from dorganism.models import Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture, OrgBatch_Image


class HiddenSimpleArrayField(forms.Field):
    widget = HiddenInput

    def clean(self, value):
        return value or []
    

class DateField(forms.Field):
    widget = forms.DateInput(attrs={'type': 'date'})
    def clean(self, value):
        return value

#=================================================================================================
# Taxonomy
#=================================================================================================
class Taxonomy_Form(forms.ModelForm):
    org_class = forms.ModelChoiceField(widget=forms.Select(attrs={'class':''}), queryset=Dictionary.objects.all(), )
    division = forms.ModelChoiceField(queryset=Dictionary.objects.all())
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['org_class'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Taxonomy.Choice_Dictionary['org_class'])]
        self.fields['division'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Taxonomy.Choice_Dictionary['division'])]

    class Meta:
        model =Taxonomy
        exclude=['urlname']
        fields=["organism_name","other_names", "code", "org_class", "tax_id", "parent_tax_id", "tax_rank", "division", "lineage" ]


class Taxonomy_Filter(Filterbase):
    organism_name = CharFilter(lookup_expr='icontains')
    lineage = CharFilter(lookup_expr='icontains')
    org_class=ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['org_class']))
    division= ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    
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

#=================================================================================================
# Organism
#=================================================================================================
class Organism_Filter(Filterbase):
    
    ID = CharFilter(field_name='organism_id', lookup_expr='icontains')
    Name = CharFilter(field_name='organism_name__organism_name', lookup_expr='icontains')
    Class = ChoiceFilter(field_name='organism_name__org_class__dict_value',  
                                      widget=forms.RadioSelect, 
                                      choices=(("GN","GN"),("GP","GP"),("FG","FG"),("MB","MB")))
    Strain = CharFilter(field_name='strain_ids', lookup_expr='icontains')
    Notes = CharFilter(field_name='strain_notes', lookup_expr='icontains')
    Type = MultipleChoiceFilter(field_name='strain_type', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    MTA = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status'], astatus__gte=0))
    Panel = MultipleChoiceFilter(field_name='strain_panel', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Type"].extra["choices"]=Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_type'], showDesc = False)
        self.filters["Panel"].extra["choices"]=Dictionary.get_aschoices(Organism.Choice_Dictionary['strain_panel'], showDesc = False)
        for i in self.filters:
            self.filters[i].label=i
   
    class Meta:
        model=Organism
        fields=[ 'Class', 'ID', 'Name','Strain',  'Notes', 'Type', 'MTA', 'Panel', ]

# -----------------------------------------------------------------
class CreateOrganism_form(forms.ModelForm):

    strain_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    # prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    oxygen_pref=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    risk_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())
    pathogen_group=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False,queryset=Dictionary.objects.all())
    mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    lab_restriction = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    res_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    gen_property=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    collect_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
   
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
       
        data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
    
        if data:
            if str(data.org_class) in ORGANISM_CLASSES:
                return data
            else:
                self.add_error('organism_name', 'Create failed with invalid organism class')
                
        else:
            self.add_error('organism_name', "Found No Organism")

        return data            


    def create_field_groups(self):
        self.group1 = [self[name] for name in Organism.FORM_GROUPS['Group1']]
        self.group2 = [self[name] for name in Organism.FORM_GROUPS['Group2']]
        self.group3 = [self[name] for name in Organism.FORM_GROUPS['Group3']]
        self.group4 = [self[name] for name in Organism.FORM_GROUPS['Group4']] 
    
    class Meta:
        model=Organism
        exclude=['organism_id',  'assoc_documents'] 

# -----------------------------------------------------------------
class UpdateOrganism_form(CreateOrganism_form):     
    class Meta:
        model=Organism
        exclude=['organism_id', 'assoc_documents'] 


   
#=================================================================================================
# Organism Batch
#=================================================================================================
class OrgBatch_Form(forms.ModelForm):

    #alphanumeric = RegexValidator(r'^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')

    # organism_id=forms.ModelChoiceField(queryset=Organism.objects.filter(astatus__gte=0), widget=forms.HiddenInput(),required=False,)
    batch_quality = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    quality_source=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    batch_id=forms.CharField(widget=forms.TextInput(attrs={'maxlength': '5', 'default':'optional input','pattern':'[0-9a-zA-Z]'}), 
                                        help_text='Optional - If empty, next number will be assigned', required=False)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)

    def __init__(self, *args, **kwargs):
        super(OrgBatch_Form, self).__init__(*args, **kwargs)
        self.fields['batch_quality'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['batch_quality'])] 
        self.fields['qc_status'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['qc_status'])] 

    class Meta:
        model =Organism_Batch
        exclude=['orgbatch_id', 'stock_level', 'organism_id']
        

class OrgBatch_UpdateForm(forms.ModelForm):

    batch_quality = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    orgbatch_id = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    quality_source=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    stock_level = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),required=False,)#SimpleArrayField(forms.IntegerField(), delimiter=';', disabled=True)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    
    def __init__(self, *args, **kwargs):   
        super().__init__(*args, **kwargs)
        instance=kwargs.get('instance')
        if instance and instance.stock_level:
            self.fields['stock_level'].initial=instance.stock_level
        self.fields['batch_quality'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['batch_quality'])]
        self.fields['qc_status'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Organism_Batch.Choice_Dictionary['qc_status'])]
        self.create_field_groups()

    def create_field_groups(self):
        self.group1 = [self[name] for name in Organism_Batch.FORM_GROUPS['Group1']]
        
    class Meta:
        model =Organism_Batch
        fields=list(model.HEADER_FIELDS.keys())
        fields+=['orgbatch_id']
        exclude=['stock_level']

# -----------------------------------------------------------------------------------    
class OrgBatch_Filter(Filterbase):
    Stock_Date = IsoDateTimeFilter(field_name='stock_date')
    class Meta:
        model = Organism_Batch
        fields = ["stock_date",  "biologist"]

#=================================================================================================
# OrgBatch Stock
#=================================================================================================

# -----------------------------------------------------------------------------------    
class OrgBatchStock_CreateForm(forms.ModelForm):

    field_order = ['orgbatch_id','stock_type', 'n_created', 'n_left', 'stock_date', 'stock_note', 'location_freezer', 'location_rack', 'location_column', 'location_slot', 'biologist']

    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))
    orgbatch_id=forms.ModelChoiceField(queryset=Organism_Batch.objects.filter(astatus__gte=0))#widget=forms.HiddenInput()
    stock_type=forms.ModelChoiceField(widget=forms.Select(attrs={'class':'', 'readonly':False}),queryset=Dictionary.objects.all(),)
    # passage_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    stock_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['stock_type'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(OrgBatch_Stock.Choice_Dictionary['stock_type'])]

    class Meta:
        model =OrgBatch_Stock
        fields='__all__'

# -----------------------------------------------------------------------------------    
class OrgBatchStock_Form(OrgBatchStock_CreateForm):
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['stock_date'].initial=self.instance.stock_date
        self.fields['n_created'].initial=self.instance.n_created

    class Meta:
        model =OrgBatch_Stock
        fields="__all__"
    
# -----------------------------------------------------------------------------------    
class OrgBatchStock_Filter(Filterbase):
    start_date = DateFilter(field_name='stock_date',lookup_expr=('gt'), widget=forms.DateInput(attrs={'type': 'date'})) 
    end_date = DateFilter(field_name='stock_date',lookup_expr=('lt'), widget=forms.DateInput(attrs={'type': 'date'}))
    Stock_Date = DateRangeFilter(field_name='stock_date')
 
    class Meta:
        model = OrgBatch_Stock
        fields = ["orgbatch_id", "Stock_Date", "start_date", "end_date"]


#=================================================================================================
# Organism Culture
#=================================================================================================
class OrgCulture_Form(forms.ModelForm):
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(widget=forms.Select(attrs={'readonly':False}), required=False,queryset=Dictionary.objects.all(),)
    culture_source= forms.ModelChoiceField(widget=forms.Select(attrs={'readonly':False}), required=False,queryset=Dictionary.objects.all(),)

    # 
    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_type'])]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_source'])]

    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys())

# -----------------------------------------------------------------------------------    
class OrgCulture_UpdateForm(forms.ModelForm):
    culture_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    culture_type= forms.ModelChoiceField(widget=forms.Select(attrs={'class':'', 'width':'fit-content','disabled': 'disabled'}), required=False,queryset=Dictionary.objects.all(),)
    culture_source= forms.ModelChoiceField(widget=forms.Select(attrs={'class':'',  'width':'fit-content','disabled': 'disabled'}), required=False,queryset=Dictionary.objects.all(),)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['culture_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_type'])]
        self.fields['culture_source'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organism_Culture.Choice_Dictionary['culture_source'])]
        self.create_field_groups()

    def create_field_groups(self):
        self.group1 = [self[name] for name in Organism_Culture.FORM_GROUPS['Group1']]

        
    class Meta:
        model =Organism_Culture
        fields=list(model.HEADER_FIELDS.keys()) 
        exclude=['culture_type', 'culture_source',]

#=================================================================================================
# OrganismBatch Image
#=================================================================================================
class OrgBatchImg_Form(forms.ModelForm):
    field_order = ['image_file','orgbatch_id','image_desc','image_source']
    image_file = forms.ImageField(label='Select an image file',required=True)
    image_type = forms.CharField(widget=forms.TextInput(attrs={'readonly':'readonly'}))
    image_name = forms.CharField(widget=forms.TextInput(attrs={'readonly':'readonly'}))
 
    def __init__(self, *args, org=None, **kwargs):
        
        super().__init__(*args, **kwargs)
        if org:
            pk = org
            organism = get_object_or_404(Organism, pk=pk)
            self.fields['orgbatch_id'].queryset = Organism_Batch.objects.filter(organism_id = organism.pk)
        
    class Meta:
        model =OrgBatch_Image
        fields="__all__"


## Stock

