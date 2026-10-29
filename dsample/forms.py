from django import forms
from applib.django.base.forms import Base_ModelForm
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
from dsample.models import  Project, ABase_Compound, COADD_Compound, ABase_Compound_Batch

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
class Project_CreateForm(Base_ModelForm):

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
        if len(Project.DETAIL_VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Project.DETAIL_VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

# -----------------------------------------------------------------
class Project_CreateMinimalForm(Base_ModelForm):

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

#=================================================================================================
# ABase Compound
#=================================================================================================
class ABase_Filter(BaseStatus_Filter):

    compound_id = CharFilter(field_name='compound_id', lookup_expr='icontains')
    compound_code = CharFilter(field_name='compound_code', lookup_expr='icontains')
    compound_name = CharFilter(field_name='compound_name', lookup_expr='icontains')
    reg_mf = CharFilter(field_name='reg_mf', lookup_expr='icontains')
    structure_type = CharFilter(field_name='structure_type', lookup_expr='icontains')
    structure_metal = CharFilter(field_name='structure_metal', lookup_expr='icontains')

    # Filter ABase_Compounds where any related ABase_Compound_Batch matches.
    chemist = ChoiceFilter(method='filter_batch_chemist', choices=[], label='Chemist')
    supplier = CharFilter(method='filter_batch_supplier', label='Supplier')
    labbook_no = CharFilter(method='filter_batch_labbook', label='LabBook No')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['chemist'].extra['choices'] = [
            (u.name, f"{u.name}") for u in ApplicationUser.objects.all().order_by('name')
        ]
        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label = self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except Exception:
                self.filters[i].label = self.filters[i].label or i

    def filter_batch_chemist(self, queryset, name, value):
        compound_ids = ABase_Compound_Batch.objects.filter(chemist_id=value).values('compound_id')
        return queryset.filter(compound_id__in=compound_ids)

    def filter_batch_supplier(self, queryset, name, value):
        compound_ids = ABase_Compound_Batch.objects.filter(supplier__icontains=value).values('compound_id')
        return queryset.filter(compound_id__in=compound_ids)

    def filter_batch_labbook(self, queryset, name, value):
        compound_ids = ABase_Compound_Batch.objects.filter(labbook_no__icontains=value).values('compound_id')
        return queryset.filter(compound_id__in=compound_ids)

    class Meta:
        model = ABase_Compound
        fields = ['compound_id', 'compound_code', 'compound_name',
                  'reg_mf', 'structure_type', 'structure_metal']

#=================================================================================================
# COADD Compound
#=================================================================================================
class COADD_Filter(BaseStatus_Filter):

    compound_id = CharFilter(field_name='compound_id', lookup_expr='icontains')
    compound_code = CharFilter(field_name='compound_code', lookup_expr='icontains')
    compound_name = CharFilter(field_name='compound_name', lookup_expr='icontains')
    std_mf = CharFilter(field_name='std_mf', lookup_expr='icontains')
    std_structure_type = CharFilter(field_name='std_structure_type', lookup_expr='icontains')
    std_metal = CharFilter(field_name='std_metal', lookup_expr='icontains')
    std_status = CharFilter(field_name='std_status', lookup_expr='icontains')
    pub_status = CharFilter(field_name='pub_status', lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label = self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except Exception:
                self.filters[i].label = i

    class Meta:
        model = COADD_Compound
        fields = ['compound_id', 'compound_code', 'compound_name',
                  'std_mf', 'std_structure_type', 'std_metal',
                  'std_status', 'pub_status']
