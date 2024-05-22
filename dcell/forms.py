from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter
from dorganism.models import  Taxonomy

from apputil.models import Dictionary, ApplicationUser, Document
from apputil.utils.filters_base import Filterbase
from adjcoadd.constants import ORGANISM_CLASSES, CELL_CLASSES #<not needed due to lack of cell classes>
from dcell.models import Cell, Cell_Batch, CellBatch_Stock



#=================================================================================================
# Cell
#=================================================================================================
class Cell_Filter(Filterbase):
    
    ID = CharFilter(field_name='cell_id', lookup_expr='icontains')
    Name = CharFilter(field_name='cell_names__cell_names', lookup_expr='icontains')
    # Class = ChoiceFilter(field_name='cell_name__org_class__dict_value',  
    #                                   widget=forms.RadioSelect, 
    #                                   choices=(("GN","GN"),("GP","GP"),("FG","FG"),("MB","MB"))) <not needed due to lack of cell classes>
    Line = CharFilter(field_name='cell_line', lookup_expr='icontains')
    Notes = CharFilter(field_name='cell_notes', lookup_expr='icontains')
    Type = MultipleChoiceFilter(field_name='cell_type', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
    MTA = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Cell.Choice_Dictionary['mta_status'], astatus__gte=0))
    Panel = MultipleChoiceFilter(field_name='cell_panel', method='multichoices_filter', 
                                             widget=forms.CheckboxSelectMultiple(attrs={'class': 'multiselect-accord'}), choices=[])
   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Type"].extra["choices"]=Dictionary.get_aschoices(Cell.Choice_Dictionary['cell_type'], showDesc = False)
        self.filters["Panel"].extra["choices"]=Dictionary.get_aschoices(Cell.Choice_Dictionary['cell_panel'], showDesc = False)
        for i in self.filters:
            self.filters[i].label=i
   
    class Meta:
        model=Cell
        fields=[ 'ID', 'Name','Line',  'Notes', 'Type', 'MTA', 'Panel', ]

# -----------------------------------------------------------------
class Cell_CreateForm(forms.ModelForm):

    cell_line= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    cell_names= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    cell_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}),required=False,)
    # prep_notes= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    mta_status = forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}),required=False, queryset=Dictionary.objects.all())
    organism_name=forms.ModelChoiceField(queryset=Taxonomy.objects.all(), widget=forms.HiddenInput(),required=False,)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    #collect_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
   
    def __init__(self, organism_name=None, *args, **kwargs): 
        self.organism_name=organism_name
        super(Cell_CreateForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['cell_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Cell.Choice_Dictionary['cell_type'], showDesc=False),)
        self.fields['cell_type'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true',})

        self.fields['cell_panel'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Cell.Choice_Dictionary['cell_panel'], showDesc=False),)
        self.fields['cell_panel'].widget.attrs.update({'class': 'form-control', 'size':'5', 'multiple': 'true'})

        self.fields['mta_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Cell.Choice_Dictionary['mta_status'])]
        self.create_field_groups()

        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                # Add the 'group-input' class to the widget attrs
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
    
    def clean_organism_name(self):       
        print(self.organism_name)
        data=get_object_or_404(Taxonomy, organism_name=self.organism_name)
    
        if data:
            if str(data.org_class) in CELL_CLASSES:
                return data
            else:
                self.add_error('org_name', 'Create failed with invalid Organism class')
                
        else:
            self.add_error('org_name', "Found No Organism")

        return data            


    def create_field_groups(self):
        self.group1 = [self[name] for name in Cell.FORM_GROUPS['Group1']]
        self.group2 = [self[name] for name in Cell.FORM_GROUPS['Group2']]
        self.group3 = [self[name] for name in Cell.FORM_GROUPS['Group3']]
        self.group4 = [self[name] for name in Cell.FORM_GROUPS['Group4']] 
    
    class Meta:
        model=Cell
        exclude=['cell_id',  'assoc_documents'] 

# -----------------------------------------------------------------
class Cell_UpdateForm(Cell_CreateForm):     
    class Meta:
        model=Cell
        exclude=['cell_id', 'assoc_documents'] 


   
#=================================================================================================
# Cell Batch
#=================================================================================================
class CellBatch_Form(forms.ModelForm):

    #alphanumeric = RegexValidator(r'^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')

    batch_id=forms.CharField(widget=forms.TextInput(attrs={'maxlength': '5', 'default':'optional input','pattern':'[0-9a-zA-Z]'}), 
                                        help_text='Optional - If empty, next number will be assigned', required=False)
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    previous_batch_id= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group'}), required=False,)
    passage_number= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group'}), required=False,)
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    batch_quality = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    quality_source=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)

    def __init__(self, *args, **kwargs):
        super(CellBatch_Form, self).__init__(*args, **kwargs)
        self.fields['batch_quality'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Cell_Batch.Choice_Dictionary['batch_quality'])] 
        self.fields['qc_status'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Cell_Batch.Choice_Dictionary['qc_status'])] 

    class Meta:
        model =Cell_Batch
        exclude=['cellbatch_id', 'stock_level', 'cell_id']
        

class CellBatch_UpdateForm(forms.ModelForm):


    cellbatch_id = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),)
    batch_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    previous_batch_id= forms.CharField(widget=forms.TextInput(), required=False,)
    passage_number= forms.CharField(widget=forms.TextInput(), required=False,)
    qc_status = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    batch_quality = forms.ModelChoiceField(required=False,queryset=Dictionary.objects.all(),)
    quality_source=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}), required=False,)
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    stock_level = forms.CharField(widget=forms.TextInput(attrs={'readonly': 'readonly'}),required=False,)
    #SimpleArrayField(forms.IntegerField(), delimiter=';', disabled=True)
    biologist=forms.ModelChoiceField(queryset=ApplicationUser.objects.all(), required=True,)
    
    def __init__(self, *args, **kwargs):   
        super().__init__(*args, **kwargs)
        instance=kwargs.get('instance')
        if instance and instance.stock_level:
            self.fields['stock_level'].initial=instance.stock_level
        self.fields['batch_quality'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Cell_Batch.Choice_Dictionary['batch_quality'])]
        self.fields['qc_status'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Cell_Batch.Choice_Dictionary['qc_status'])]
        self.create_field_groups()

    def create_field_groups(self):
        self.group1 = [self[name] for name in Cell_Batch.FORM_GROUPS['Group1']]
        
    class Meta:
        model =Cell_Batch
        fields=list(model.HEADER_FIELDS.keys())
        fields+=['cellbatch_id']
        exclude=['stock_level']

# -----------------------------------------------------------------------------------    
class CellBatch_Filter(Filterbase):
    Stock_Date = IsoDateTimeFilter(field_name='stock_date')
    class Meta:
        model = Cell_Batch
        fields = ["stock_date",  "biologist"]

#=================================================================================================
# CellBatch Stock
#=================================================================================================

# -----------------------------------------------------------------------------------    
class CellBatchStock_CreateForm(forms.ModelForm):

    field_order = ['cellbatch_id','stock_type', 'n_created', 'n_left', 'stock_date', 'stock_note', 'location_freezer', 'location_rack', 'location_column', 'location_slot', 'biologist']

    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))
    cellbatch_id=forms.ModelChoiceField(queryset=Cell_Batch.objects.filter(astatus__gte=0))#widget=forms.HiddenInput()
    stock_type=forms.ModelChoiceField(widget=forms.Select(attrs={'class':'', 'readonly':False}),queryset=Dictionary.objects.all(),)
    # passage_notes=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)
    stock_note=forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '3'}), required=False,)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['stock_type'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(CellBatch_Stock.Choice_Dictionary['stock_type'])]

    class Meta:
        model =CellBatch_Stock
        fields='__all__'

# -----------------------------------------------------------------------------------    
class CellBatchStock_Form(CellBatchStock_CreateForm):
    stock_date=forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}))
    n_created=forms.IntegerField(widget=forms.NumberInput(attrs={'type': 'number'}))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['stock_date'].initial=self.instance.stock_date
        self.fields['n_created'].initial=self.instance.n_created

    class Meta:
        model =CellBatch_Stock
        fields="__all__"
    
# -----------------------------------------------------------------------------------    
class CellBatchStock_Filter(Filterbase):
    start_date = DateFilter(field_name='stock_date',lookup_expr=('gt'), widget=forms.DateInput(attrs={'type': 'date'})) 
    end_date = DateFilter(field_name='stock_date',lookup_expr=('lt'), widget=forms.DateInput(attrs={'type': 'date'}))
    Stock_Date = DateRangeFilter(field_name='stock_date')
 
    class Meta:
        model = CellBatch_Stock
        fields = ["cellbatch_id", "Stock_Date", "start_date", "end_date"]
