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
from applib.django.filters import BaseStatus_Filter
 

#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50
#from dsummary.models import Summary_ScreenRun

#=================================================================================================
# Screen_Run
#=================================================================================================
# class SumScreenRun_Filter(BaseStatus_Filter):
    
#     f_RunID = CharFilter(field_name='run_id__run_id', lookup_expr='icontains', label="Run ID")
#     f_RunName = CharFilter(field_name='run_id__run_name', lookup_expr='icontains', label="Run Name")
#     f_RunType=ChoiceFilter(field_name='run_id__run_type',widget=forms.RadioSelect, choices=[], label="Run Type")
#     f_RunStatus=ChoiceFilter(field_name='run_id__run_status',widget=forms.RadioSelect, choices=[], label="Run Status")
    
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.filters["f_RunType"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
#         self.filters["f_RunStatus"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

#         # Set Filter label to the Fields VerboseName or Filter Name
#         # for i in self.filters:
#         #     try:
#         #         self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
#         #     except:
#         #         self.filters[i].label=i
#     class Meta:
#         model=Summary_ScreenRun
#         fields=[ 'f_RunID', 'f_RunName','f_RunType','f_RunStatus']

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
        
    class Meta:
        model=Screen_Run
        exclude=Screen_Run.CALCULATED_FIELDS

    def create_field_groups(self):
        if len(Screen_Run.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Screen_Run.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

# -----------------------------------------------------------------
class ScreenRun_UpdateForm(ScreenRun_CreateForm):     
    class Meta:
        model=Screen_Run
        exclude=['run_id']+Screen_Run.CALCULATED_FIELDS
