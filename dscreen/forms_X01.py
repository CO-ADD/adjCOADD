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
from apputil.utils.filters_base import Filterbase
#from adjcoadd.constants import ORGANISM_CLASSES


#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50

#=================================================================================================
# Screen_Run
#=================================================================================================
class ScreenRun_Filter(Filterbase):
    
    Run_Type=ChoiceFilter(field_name='run_type',widget=forms.RadioSelect, choices=[], empty_label=None)
    Run_Status=ChoiceFilter(field_name='run_status',widget=forms.RadioSelect, choices=[], empty_label=None)
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["Run_Type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        self.filters["Run_Status"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]
        # for i in self.filters:
        #     self.filters[i].label=i
   
    class Meta:
        model=Screen_Run
        fields=[ 'run_id', 'run_name','Run_Type','Run_Status']

# -----------------------------------------------------------------
class ScreenRun_Form(forms.ModelForm):

    run_type = forms.ChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, choices= [], )
    run_status = forms.ChoiceField(widget=forms.Select(attrs={'class':'form-select'}), required=False, choices= [], )
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        self.fields['run_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])]
        self.fields['run_status'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_status'])]

        # Create groups of fields for View 
        self.create_field_groups()

        # Add the 'group-input' class to the widget attrs
        for field in self.fields.values():
            if isinstance(field.widget, forms.TextInput) or isinstance(field.widget, forms.NumberInput):
                attrs = field.widget.attrs
                attrs['class'] = attrs.get('class', '') + 'input-group'
                field.widget.attrs = attrs
                
    def create_field_groups(self):
        if len(Screen_Run.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Screen_Run.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])
    class Meta:
        model =Screen_Run
        fields='__all__'
        exclude = []
        
    
# -----------------------------------------------------------------
class ScreenRun_CreateForm(forms.ModelForm):
    run_id=forms.CharField(widget=forms.TextInput(attrs={'maxlength': '15', 'default':'optional input','pattern':'[0-9a-zA-Z]'}), 
                                        help_text='Optional - If empty, next ID will be assigned as per RunType', required=False)

    def __init__(self, *args, **kwargs):
        super(ScreenRun_CreateForm, self).__init__(*args, **kwargs)
        self.fields['run_type'].choices=[(obj.dict_value, repr(obj)) for obj in Dictionary.get_filterobj(Screen_Run.DICTIONARY_FIELDS['run_type'])] 
        self.fields['run_status'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Screen_Run.DICTIONARY_FIELDS['run_status'], showDesc=False),)

    class Meta:
        model=Screen_Run
        exclude=Screen_Run.CALCULATED_FIELDS 
        
# -----------------------------------------------------------------
class ScreenRun_DetailForm(forms.ModelForm):

    def __init__(self, organism_name=None, *args, **kwargs): 
        super(ScreenRun_DetailForm, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name
        self.fields['run_type'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Screen_Run.DICTIONARY_FIELDS['run_type'], showDesc=False),)
        self.fields['run_status'].widget = forms.SelectMultiple(choices = Dictionary.get_aschoices(Screen_Run.DICTIONARY_FIELDS['run_status'], showDesc=False),)

        self.helper = FormHelper(self)
        self.helper.form_id = 'ScreenUpdateForm'


    class Meta:
        model=Screen_Run
        exclude=Screen_Run.CALCULATED_FIELDS 

# -----------------------------------------------------------------
class ScreenRun_UpdateForm(ScreenRun_DetailForm):     
    class Meta:
        model=Screen_Run
        exclude=['run_id']  
