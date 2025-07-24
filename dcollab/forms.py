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
from django_countries.data import COUNTRIES 

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
    Country = ChoiceFilter(field_name='country', choices=CountryField().choices)
    #Country = ChoiceFilter(field_name='country', choices=sorted(COUNTRIES.items()))
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters["organisation_type"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]
        #self.filters['Country'].extra["choices"] = self.Meta.model.get_field_choices(field_name='country__name')

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
    organisation_code= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    organisation_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    country = CountryField()
    organisation_type=ChoiceFilter(field_name='organisation_type',choices=[], empty_label=None)

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

# -----------------------------------------------------------------
class Organisation_UpdateForm(Organisation_CreateForm): 
    organisation_id = forms.CharField(disabled=True)  
    
    class Meta:
        model=Organisation
        exclude=[]

    field_order = Organisation.LIST_VIEW_FIELDS

#=================================================================================================
# Collab User
#=================================================================================================
class CollabUser_Filter(BaseStatus_Filter):

    FILTERSET_DICT = {
        'Organisation':      {'lookup':'choice','field_name':'organisation_id__organisation_name'},
    }


    Organisation = ChoiceFilter(field_name='organisation_id__organisation_name', choices=[], label="Organisation")
    Country = ChoiceFilter(field_name='country', choices=CountryField().choices)
    #Country = ChoiceFilter(field_name='country', choices=sorted(COUNTRIES.items()))
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Extract Filter Dictionary
        _filter_dict = {}
        if 'filterset_dict' in kwargs:
            print(kwargs['filterset_dict'])
            for _key, _item in self.FILTERSET_DICT.items():
                if _key in kwargs['filterset_dict']:
                    _filter_dict[_item['field_name']] = kwargs['filterset_dict'][_key][0]
            kwargs.pop('filterset_dict')
            
        # Initialise FilterSet and Choices
        super().__init__(*args, **kwargs)
        for _key, _item in self.FILTERSET_DICT.items():
            if _item['lookup'] == 'choice':
                self.filters[_key].extra["choices"] = self.Meta.model.get_field_choices(field_name=_item['field_name'],filter_dict=_filter_dict)

        #self.filters["Organisation"].extra['choices']=[(obj.dict_value, str(obj)) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]
        #self.filters['Country'].extra["choices"] = self.Meta.model.get_field_choices(field_name='country__name')
        #self.filters['Organisation'].extra["choices"] = self.Meta.model.get_field_choices(field_name='organisation_id__organisation_name',filter_dict=_filter_dict)

        # Set Filter label to the Fields VerboseName or Filter Name
        # for i in self.filters:
        #     try:
        #         self.filters[i].label=self.Meta.model._meta.get_field(self.filters[i].field_name).verbose_name
        #     except:
        #         self.filters[i].label=i
    
    class Meta:
        model=Collab_User
        fields=['first_name','last_name','email','Organisation','Country',]

    # user_id = models.CharField(max_length=15, primary_key=True, verbose_name = "User ID")
    # title = models.CharField(max_length=15, blank=True, verbose_name = "Title")
    # first_name = models.CharField(max_length=50, blank=True, verbose_name = "First Code")
    # last_name = models.CharField(max_length=50, blank=True, verbose_name = "Last Code")
    # position = models.CharField(max_length=100, blank=True, verbose_name = "Position")

    # email1 = models.EmailField(max_length=254, blank=True, verbose_name = "EMail 1")
    # email2 = models.EmailField(max_length=254, blank=True, verbose_name = "EMail 2")
    # active_email = models.SmallIntegerField(default=0, blank=True, verbose_name ="Active")

    # phone = models.CharField(max_length=50, blank=True, verbose_name = "Phone")
    # subscribed = models.BooleanField(default=False, blank=True, verbose_name = "Newsletter")
    # portal_userid = models.CharField(max_length=50, blank=True, verbose_name = "Portal UserID")
    # portal_pw = models.CharField(max_length=50, blank=True, verbose_name = "Portal Password")

    # organisation_id = models.ForeignKey(Organisation, null=True, blank=True, verbose_name = "Organisation ID", on_delete=models.DO_NOTHING,
    #     db_column="organisation_id", related_name="%(class)s_organisation_id")
        
    # department = models.CharField(max_length=250, blank=True, verbose_name = "Department")
    # postal_address = models.CharField(max_length=250, blank=True, verbose_name = "Postal Address")
    # city = models.CharField(max_length=250, blank=True, verbose_name = "City")
    # country = CountryField(verbose_name = "Country")



# -----------------------------------------------------------------
class CollabUser_CreateForm(forms.ModelForm):

    # PK to add help text
    #organisation_code= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    #organisation_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    organisation = ChoiceFilter(field_name='organisation_id__organisation_name', choices=[], label="Organisation")
    country = CountryField()
    #organisation_type=ChoiceFilter(field_name='organisation_type',choices=[], empty_label=None)

    def __init__(self, *args, **kwargs): 
        super(CollabUser_CreateForm, self).__init__(*args, **kwargs)
        
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        #self.fields['organisation_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]

        self.create_field_groups()

    class Meta:
        model=Collab_User
        exclude=['user_id']

    def create_field_groups(self):
        if len(Collab_User.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Collab_User.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   

#=================================================================================================
# Collab Group
#=================================================================================================
class CollabGroup_Filter(BaseStatus_Filter):
    
    FILTERSET_DICT = {
        'Organisation':      {'lookup':'choice','field_name':'organisation_id__organisation_name'},
    }


    Organisation = ChoiceFilter(field_name='organisation_id__organisation_name', choices=[], label="Organisation")
    Country = ChoiceFilter(field_name='country', choices=CountryField().choices)
    mta_status = ModelChoiceFilter(field_name='mta_status', queryset=Dictionary.objects.filter(dict_class=Collab_Group.DICTIONARY_FIELDS['mta_status'], astatus__gte=0))
    #Country = ChoiceFilter(field_name='country', choices=sorted(COUNTRIES.items()))
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Extract Filter Dictionary
        _filter_dict = {}
        if 'filterset_dict' in kwargs:
            print(kwargs['filterset_dict'])
            for _key, _item in self.FILTERSET_DICT.items():
                if _key in kwargs['filterset_dict']:
                    _filter_dict[_item['field_name']] = kwargs['filterset_dict'][_key][0]
            kwargs.pop('filterset_dict')
            
        # Initialise FilterSet and Choices
        super().__init__(*args, **kwargs)
        for _key, _item in self.FILTERSET_DICT.items():
            if _item['lookup'] == 'choice':
                self.filters[_key].extra["choices"] = self.Meta.model.get_field_choices(field_name=_item['field_name'],filter_dict=_filter_dict)

    class Meta:
        model=Collab_Group
        fields=['group_code','Organisation','Country','mta_status']

# -----------------------------------------------------------------
class CollabGroup_CreateForm(forms.ModelForm):

    # PK to add help text
    organisation_code= forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    organisation_name = forms.CharField(widget=forms.Textarea(attrs={'class': 'input-group', 'rows': '2'}),required=False,)
    country = CountryField()
    organisation_type=ChoiceFilter(field_name='organisation_type',choices=[], empty_label=None)

    def __init__(self, *args, **kwargs): 
        super(CollabGroup_CreateForm, self).__init__(*args, **kwargs)
        
        # Set Labels from Model Definitions
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

        # Set Dictionary values
        self.fields['organisation_type'].choices=[(obj.dict_value, obj.strtml()) for obj in Dictionary.get_filterobj(Organisation.DICTIONARY_FIELDS['organisation_type'])]

        self.create_field_groups()

    class Meta:
        model=Collab_Group
        exclude=['group_id']

    def create_field_groups(self):
        if len(Collab_Group.VIEW_GROUPS) > 0:
            self.groups = []
            for grp in Collab_Group.VIEW_GROUPS:
                self.groups.append([self[name] for name in grp])   


