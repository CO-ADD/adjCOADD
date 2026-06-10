from django import forms
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator
from django.forms import ModelForm
from django.shortcuts import get_object_or_404
from django.forms.widgets import HiddenInput
from django.contrib.postgres.forms import SimpleArrayField, SplitArrayField
from django_filters import DateRangeFilter, CharFilter, ModelChoiceFilter, ChoiceFilter, MultipleChoiceFilter, IsoDateTimeFilter, DateFromToRangeFilter, DateFilter

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Fieldset, Submit
from django_countries.fields import CountryField

from apputil.models import Dictionary, ApplicationUser, Document
from adjcoadd.constants import PROJECT_COMPOUND_STATUS, PROJECT_SCREEN_STATUS, PROJECT_DATA_STATUS, PROJECT_REPORT_STATUS
from applib.django.base.filters import BaseStatus_Filter, TrigramFilter
 
#dSample
from dsample.models import  Project

#=================================================================================================
# Project
#=================================================================================================
class Project_Filter(BaseStatus_Filter):
    
    project_type=ChoiceFilter(field_name='project_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    project_status=ChoiceFilter(field_name='project_status',widget=forms.RadioSelect, choices=[], empty_label=None)
    Country = ChoiceFilter(field_name='group_id__country', choices=CountryField().choices,)
    #Organisation = CharFilter(field_name='group_id__organisation_id__organisation_name', lookup_expr='icontains',label='Organisation')
    Organisation = TrigramFilter(field_name='group_id__organisation_id__organisation_name', label='Organisation', byword=True)
    
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["project_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_type'])]
        self.filters["project_status"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_status'])]

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
                
    class Meta:
        model=Project
        fields=[ 'project_id','project_name','project_type','project_status','Country','Organisation']
        exclude = ['group_id.organisation_id.organisation_name',
                   ]


# -----------------------------------------------------------------
class Project_CreateForm(forms.ModelForm):

    # PK to add help text
    project_id= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}),required=False,help_text="Leave empty to use next Pnnnnn number")

    # DateFields
    pub_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
    received = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
    completed = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)

    compound_status = forms.MultipleChoiceField(choices=[(i,i) for i in PROJECT_COMPOUND_STATUS], widget=forms.CheckboxSelectMultiple(),required=False)
    screen_status = forms.MultipleChoiceField(choices=[(i,i) for i in PROJECT_SCREEN_STATUS], widget=forms.CheckboxSelectMultiple(),required=False)
    data_status = forms.MultipleChoiceField(choices=[(i,i) for i in PROJECT_DATA_STATUS], widget=forms.CheckboxSelectMultiple(),required=False)
    report_status = forms.MultipleChoiceField(choices=[(i,i) for i in PROJECT_REPORT_STATUS], widget=forms.CheckboxSelectMultiple(),required=False)
 
    # Simple Array Fields
    stock_status = SimpleArrayField(forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'})), required=False, delimiter=';', max_length=20)
 
    # TextFields - 2 rows (Normal,short CharFields do not need definition)
    project_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    provided_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    project_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    stock_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    compound_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    screen_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    data_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    report_comment= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    source= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    reference= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)

    def __init__(self, *args, **kwargs): 
        super(Project_CreateForm, self).__init__(*args, **kwargs)
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        self.fields['project_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_type'])]
        self.fields['project_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_status'])]
        self.fields['pub_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['pub_status'])]
        self.fields['provided_container'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['provided_container'])]
        self.fields['stock_conc_unit'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['stock_conc_unit'])]

        # Make Calculated fields ReadOnly and Hidden
        for field_name in Project.CALCULATED_FIELDS:
            self.fields[field_name].widget.attrs['readonly'] = 'readonly'
            self.fields[field_name].widget = self.fields[field_name].hidden_widget()

        self.create_field_groups()

    class Meta:
        model=Project
        exclude=Project.ORACLE_FIELDS

    def create_field_groups(self):
        if len(Project.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Project.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

# -----------------------------------------------------------------
class Project_CreateMinimalForm(forms.ModelForm):

    # PK to add help text
    project_id= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}),required=False,help_text="Leave empty to use next Pnnnnn number")
    #project_name= forms.CharField(widget=forms.TextInput(attrs={'class': 'input-group'}),required=False,initial="Project Name")
    # DateFields
    #pub_date = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)
    received = forms.DateField(widget=forms.DateInput(attrs={'type': 'date'}), required=False)

    def __init__(self, *args, **kwargs): 
        super(Project_CreateMinimalForm, self).__init__(*args, **kwargs)

       # Set Dictionary values
        self.fields['project_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_type'])]
        self.fields['project_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['project_status'])]
        self.fields['pub_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['pub_status'])]
        self.fields['provided_container'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['provided_container'])]
        #self.fields['stock_conc_unit'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Project.DICTIONARY_FIELDS['stock_conc_unit'])]

    class Meta:
        model=Project
        fields=['project_id','project_name','project_comment','project_type','project_status', 'received',
                'pub_status','provided_container',
                #'group_id',
                #'project_members',
                ]
    

# -----------------------------------------------------------------
class Project_UpdateForm(Project_CreateForm):     
    class Meta:
        model=Project
        exclude=['project_id','project_members'] + Project.ORACLE_FIELDS 
