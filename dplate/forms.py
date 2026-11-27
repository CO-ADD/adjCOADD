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
from applib.django.base.filters import BaseStatus_Filter
 
#dSample
from dplate.models import  TestPlate, MasterPlate, Labware

#=================================================================================================
# TestPlate
#=================================================================================================
class TestPlate_Filter(BaseStatus_Filter):

    plate_id = CharFilter(field_name='plate_id', lookup_expr='icontains', label='Plate ID')
    run_id = CharFilter(field_name='run_id__run_id', lookup_expr='icontains', label='Run ID')
    assay_id = CharFilter(field_name='assay_id__assay_id', lookup_expr='icontains', label='Assay ID')
    result_type = ChoiceFilter(field_name='result_type', widget=forms.RadioSelect, choices=[], empty_label=None)
    plate_quality = ChoiceFilter(field_name='plate_quality', widget=forms.RadioSelect, choices=[], empty_label=None)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["result_type"].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['result_type'])]
        self.filters["plate_quality"].extra['choices'] = [(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['plate_quality'])]

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label = self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except Exception:
                self.filters[i].label = self.filters[i].label or i

    class Meta:
        model = TestPlate
        fields = ['plate_id', 'run_id', 'result_type', 'assay_id', 'plate_quality']


#=================================================================================================
# MasterPlate
#=================================================================================================
class MasterPlate_Filter(BaseStatus_Filter):

    plate_id = CharFilter(field_name='plate_id', lookup_expr='icontains', label='Plate ID')
    run_id = CharFilter(field_name='run_id__run_id', lookup_expr='icontains', label='Run ID')
    plate_quality = ChoiceFilter(field_name='plate_quality', widget=forms.RadioSelect, choices=[], empty_label=None)
    prep_date = DateFilter(field_name='prep_date', widget=forms.DateInput(attrs={'type': 'date'}), label='Prep Date')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['plate_quality'].extra['choices'] = [
            (obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(MasterPlate.DICTIONARY_FIELDS['plate_quality'])
        ]
        for i in self.filters:
            try:
                self.filters[i].label = self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except Exception:
                self.filters[i].label = self.filters[i].label or i

    class Meta:
        model = MasterPlate
        fields = ['plate_id', 'run_id', 'plate_quality', 'prep_date']


#=================================================================================================
# Labware
#=================================================================================================
class Labware_Filter(BaseStatus_Filter):

    labware_id = CharFilter(field_name='labware_id', lookup_expr='icontains', label='Labware ID')
    labware_name = CharFilter(field_name='labware_name', lookup_expr='icontains', label='Name')
    labware_type = CharFilter(field_name='labware_type', lookup_expr='icontains', label='Type')
    brand = CharFilter(field_name='brand', lookup_expr='icontains', label='Brand')

    class Meta:
        model = Labware
        fields = ['labware_id', 'labware_name', 'labware_type', 'brand']


#=================================================================================================
# Update Forms
#=================================================================================================
class TestPlate_UpdateForm(Base_ModelForm):
    plate_id = forms.CharField(disabled=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['plate_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['plate_type'])]
        self.fields['result_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['result_type'])]
        self.fields['plate_quality'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['plate_quality'])]

    class Meta:
        model = TestPlate
        fields = ['plate_id', 'plate_type', 'labware_id', 'run_id', 'assay_id', 'result_type',
                  'plate_quality', 'test_date', 'control_layout', 'plate_comment', 'zfactor', 'process_status']
        field_order = ['plate_id', 'plate_type', 'labware_id', 'run_id', 'assay_id', 'result_type',
                       'plate_quality', 'test_date', 'control_layout', 'plate_comment', 'zfactor', 'process_status']


class MasterPlate_UpdateForm(Base_ModelForm):
    plate_id = forms.CharField(disabled=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['plate_type'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(MasterPlate.DICTIONARY_FIELDS['plate_type'])]
        self.fields['plate_quality'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(MasterPlate.DICTIONARY_FIELDS['plate_quality'])]

    class Meta:
        model = MasterPlate
        fields = ['plate_id', 'plate_type', 'labware_id', 'run_id', 'plate_quality', 'prep_date',
                  'plating', 'dilution_layout', 'master_solvent', 'master_conc', 'master_volume', 'plate_comment', 'process_status']
        field_order = ['plate_id', 'plate_type', 'labware_id', 'run_id', 'plate_quality', 'prep_date',
                       'plating', 'dilution_layout', 'master_solvent', 'master_conc', 'master_volume', 'plate_comment', 'process_status']


class Labware_Form(Base_ModelForm):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['plate_material'].choices = [(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Labware.DICTIONARY_FIELDS['plate_material'])]

    class Meta:
        model = Labware
        fields = ['labware_id', 'labware_name', 'labware_type', 'plate_size', 'plate_type', 'plate_color',
                  'plate_material', 'well_shape', 'well_bottom', 'well_type', 'well_size',
                  'working_volume', 'brand', 'model']
        field_order = ['labware_id', 'labware_name', 'labware_type', 'plate_size', 'plate_type', 'plate_color',
                       'plate_material', 'well_shape', 'well_bottom', 'well_type', 'well_size',
                       'working_volume', 'brand', 'model']


class Labware_UpdateForm(Labware_Form):
    labware_id = forms.CharField(disabled=True)

    class Meta:
        model = Labware
        fields = ['labware_id', 'labware_name', 'labware_type', 'plate_size', 'plate_type', 'plate_color',
                  'plate_material', 'well_shape', 'well_bottom', 'well_type', 'well_size',
                  'working_volume', 'brand', 'model']
        field_order = ['labware_id', 'labware_name', 'labware_type', 'plate_size', 'plate_type', 'plate_color',
                       'plate_material', 'well_shape', 'well_bottom', 'well_type', 'well_size',
                       'working_volume', 'brand', 'model']
