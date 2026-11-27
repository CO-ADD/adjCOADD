from django import forms
from applib.django.base.forms import Base_ModelForm
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Fieldset, Submit

from apputil.models import Dictionary, ApplicationUser, Document
from applib.django.base.filters import BaseStatus_Filter
 

#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50
#from dsummary.models import Summary_ScreenRun

#=================================================================================================
# Screen_Run
#=================================================================================================
class ScreenRun_Filter(BaseStatus_Filter):
    
    run_type=ChoiceFilter(field_name='run_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    run_status=ChoiceFilter(field_name='run_status',widget=forms.RadioSelect, choices=[], empty_label=None)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["run_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        self.filters["run_status"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=Screen_Run
        fields=[ 'run_id', 'run_name','run_type','run_status']

# -----------------------------------------------------------------
class ScreenRun_CreateForm(Base_ModelForm):

    # PK to add help text
    run_id= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}),required=False,help_text="Leave empty to use next PSR/HCR/.. number")
    #run_type=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=True,queryset=Dictionary.objects.all())
    #run_status=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())

    # DateFields
    run_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)

    # TextFields - 2 rows (Normal,short CharFields do not need definition)
    assay_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    run_conditions= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    run_issues= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    run_project= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)

    def __init__(self, *args, **kwargs): 
        super().__init__(*args, **kwargs)
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        self.fields['run_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        self.fields['run_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

        # Additional attributes
        # self.fields["run_id"].widget.attrs.update({"class":"special"})
        # self.fields["run_id"].widget.attrs.update(size=40)

        # Create groups of fields for View 
        self.create_field_groups()
        
        # Add the 'group-input' class to the widget attrs
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
            
        # Make Calculated fields ReadOnly and Hidden
        for field_name in Screen_Run.CALCULATED_FIELDS:
            self.fields[field_name].widget.attrs['readonly'] = True
            self.fields[field_name].widget = self.fields[field_name].hidden_widget()
            
    class Meta:
        model=Screen_Run
        exclude=[]
        #exclude=Screen_Run.CALCULATED_FIELDS

    def create_field_groups(self):
        if len(Screen_Run.DETAIL_VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Screen_Run.DETAIL_VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

# -----------------------------------------------------------------
class ScreenRun_UpdateForm(ScreenRun_CreateForm):

    def __init__(self, *args, **kwargs): 
        super(ScreenRun_UpdateForm, self).__init__(*args, **kwargs)

        # Make Calculated fields ReadOnly
        for field in Screen_Run.CALCULATED_FIELDS:
            self.fields[field].widget.attrs['readonly'] = True
    class Meta:
        model=Screen_Run
        exclude=['run_id']


#=================================================================================================
# Assay
#=================================================================================================
class Assay_Filter(BaseStatus_Filter):
    
    # run_type=ChoiceFilter(field_name='run_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    # run_status=ChoiceFilter(field_name='run_status',widget=forms.RadioSelect, choices=[], empty_label=None)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self.filters["run_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        # self.filters["run_status"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=Assay
        fields=[ 'assay_id', 'assay_subtype','assay_code']

#----------------------------------------------------------------
class Assay_CreateForm(Base_ModelForm):

    def __init__(self, *args, **kwargs): 
        super().__init__(*args, **kwargs)
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name


        # Create groups of fields for View 
        self.create_field_groups()
        
        # Add the 'group-input' class to the widget attrs
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
        
        # Make Calculated fields ReadOnly and Change to Hidden Widget
        # for field in Assay.CALCULATED_FIELDS:
        #     self.fields[field].widget.attrs['readonly'] = True
        #     self.fields[field].widget = self.fields[field].hidden_widget()
                    
    class Meta:
        model=Assay
        exclude=['assay_id']
        #exclude=Screen_Run.CALCULATED_FIELDS

    def create_field_groups(self):
        if len(Assay.DETAIL_VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Assay.DETAIL_VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

# -----------------------------------------------------------------
class Assay_UpdateForm(Assay_CreateForm):

    def __init__(self, *args, **kwargs): 
        super(Assay_UpdateForm, self).__init__(*args, **kwargs)

        # Make Calculated fields ReadOnly
        # for field in Assay.CALCULATED_FIELDS:
        #     self.fields[field].widget.attrs['readonly'] = True
    class Meta:
        model=Assay
        exclude=['assay_id']


#=================================================================================================
# AssayData (MIC / CC50 / HC50)
#=================================================================================================
class AssayData_BaseFilter(BaseStatus_Filter):
    """
    Shared filter fields for the (very large) AssayData result tables.
    """

    testplate_id = CharFilter(field_name='testplate_id__plate_id', lookup_expr='icontains', label='Plate')
    testwell_id = CharFilter(field_name='testwell_id', lookup_expr='icontains', label='Well')
    run_id = CharFilter(field_name='run_id__run_id', lookup_expr='icontains', label='Run ID')
    assay_id = CharFilter(field_name='assay_id__assay_id', lookup_expr='icontains', label='Assay ID')
    cmpbatch_id = CharFilter(field_name='cmpbatch_id__cmpbatch_id', lookup_expr='icontains', label='CmpBatch')
    act_type = CharFilter(field_name='act_type', lookup_expr='icontains', label='Act Type')
    data_quality = CharFilter(field_name='data_quality__dict_value', lookup_expr='icontains', label='Quality')

    class Meta:
        model = None
        fields = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Keep the explicit labels above (FK paths cannot be resolved by
        # get_field so do not attempt to overwrite them).


class AssayData_MIC_Filter(AssayData_BaseFilter):

    class Meta:
        model = AssayData_MIC
        fields = ['testplate_id', 'testwell_id', 'run_id', 'assay_id',
                  'cmpbatch_id', 'act_type', 'data_quality']


class AssayData_CC50_Filter(AssayData_BaseFilter):

    class Meta:
        model = AssayData_CC50
        fields = ['testplate_id', 'testwell_id', 'run_id', 'assay_id',
                  'cmpbatch_id', 'act_type', 'data_quality']


class AssayData_HC50_Filter(AssayData_BaseFilter):

    class Meta:
        model = AssayData_HC50
        fields = ['testplate_id', 'testwell_id', 'run_id', 'assay_id',
                  'cmpbatch_id', 'act_type', 'data_quality']