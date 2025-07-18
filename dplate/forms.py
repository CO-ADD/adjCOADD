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
from applib.django.base.filters import BaseStatus_Filter
 
#dSample
from dplate.models import  TestPlate

#=================================================================================================
# TestPlate
#=================================================================================================
class TestPlate_Filter(BaseStatus_Filter):
    
    result_type=ChoiceFilter(field_name='result_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    plate_quality=ChoiceFilter(field_name='plate_quality',widget=forms.RadioSelect, choices=[], empty_label=None)
    run_id = ChoiceFilter(field_name='run_id', widget=forms.RadioSelect, choices=[], empty_label=None)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["result_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['result_type'])]
        self.filters["plate_quality"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['plate_quality'])]
        self.filters['run_id'].extra["choices"] = self.Meta.model.get_field_choices(field_name='run_id')

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=TestPlate
        fields=[ 'plate_id','run_id','result_type','assay_id','plate_quality',]
