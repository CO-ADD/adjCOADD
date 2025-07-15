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
 
#-- dCollab --------------------------------------------------------------------
from dcollab.models import Organisation, Collab_Group, Collab_User

#=================================================================================================
# Organisation
#=================================================================================================
class Organisation_Filter(BaseStatus_Filter):
    
    organisation_type=ChoiceFilter(field_name='organisation_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Country = ChoiceFilter(field_name='organisation_id__country', choices=CountryField().choices,)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["organisation_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]
        #self.filters['Country'].extra["choices"] = self.Meta.model.get_field_choices(field_name='group_id__country')

        # Set Filter label to the Fields VerboseName or Filter Name
        for i in self.filters:
            try:
                self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
            except:
                self.filters[i].label=i
    class Meta:
        model=Organisation
        fields=[ 'organisation_id','organisation_name','organisation_type','organisation_code','Country',]

# -----------------------------------------------------------------
class Organisation_CreateForm(forms.ModelForm):

    # PK to add help text
    organisation_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    organisation_code= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    country = CountryField()

    def __init__(self, *args, **kwargs): 
        super(Organisation_CreateForm, self).__init__(*args, **kwargs)
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        self.fields['organisation_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]

        self.create_field_groups()

    class Meta:
        model=Organisation
        exclude=['organisation_id']

    def create_field_groups(self):
        if len(Organisation.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Organisation.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   
