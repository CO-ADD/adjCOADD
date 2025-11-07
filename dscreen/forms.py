from django import forms
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
class ScreenRun_CreateForm(forms.ModelForm):

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
        super(ScreenRun_CreateForm, self).__init__(*args, **kwargs)
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
        
        # Make Calculated fields ReadOnly
        for field in Screen_Run.CALCULATED_FIELDS:
            self.fields[field].widget.attrs['readonly'] = True

        
    class Meta:
        model=Screen_Run
        exclude=[]
        #exclude=Screen_Run.CALCULATED_FIELDS

    def create_field_groups(self):
        if len(Screen_Run.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Screen_Run.VIEW_GROUPS:
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

class Assay_CreateForm(forms.ModelForm):

    # PK to add help text
    #run_id= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}),required=False,help_text="Leave empty to use next PSR/HCR/.. number")
    #run_type=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=True,queryset=Dictionary.objects.all())
    #run_status=forms.ModelChoiceField(widget=forms.Select(attrs={'class': 'form-control'}), required=False,queryset=Dictionary.objects.all())

    # DateFields
    #run_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)

    # TextFields - 2 rows (Normal,short CharFields do not need definition)
    # assay_note= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    # test_media= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    # run_issues= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    # run_project= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)

    def __init__(self, *args, **kwargs): 
        super().__init__(*args, **kwargs)
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        # self.fields['run_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        # self.fields['run_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

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
        
        # Make Calculated fields ReadOnly
        # for field in Assay.CALCULATED_FIELDS:
        #     self.fields[field].widget.attrs['readonly'] = True

        
    class Meta:
        model=Assay
        exclude=['assay_id']
        #exclude=Screen_Run.CALCULATED_FIELDS

    def create_field_groups(self):
        if len(Assay.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Assay.VIEW_GROUPS:
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