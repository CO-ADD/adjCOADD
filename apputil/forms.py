#import django_filters
from django import forms
from .models import  ApplicationUser, ApplicationLog, Dictionary, Document
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate
from django.shortcuts import get_object_or_404
from django_filters import DateRangeFilter, CharFilter, BooleanFilter, ChoiceFilter, DateFromToRangeFilter, DateFilter

from apputil.utils.filters_base import Filterbase, Filterbase_base


#------------------------------------------------------------------------
# --Login Form--
class Login_Form(AuthenticationForm):
#------------------------------------------------------------------------
    username=forms.CharField(widget=forms.TextInput(attrs={'class':'form-control',  'id': 'user-input', 'autocomplete':'off'}), label='user-input')
    password= forms.CharField(widget=forms.PasswordInput(
    attrs={'class':'form-control', 'id': 'password-input','type':'password', 'name': 'password', 'autocomplete':'new-password'}),
    label='password-input')
    def clean(self):
        username = self.cleaned_data.get('username')
        password = self.cleaned_data.get('password')
        if ApplicationUser.objects.filter(username=username):
            self.user_cache = authenticate(username=username,
										   password=password)
        else:
            self.user_cache=None
        if self.user_cache is None:
            raise forms.ValidationError(
					self.error_messages['invalid_login'],
					code='user not existed',
					params={'username': self.username_field.verbose_name},
				)
        # for someone deleted to not an application user 
        elif get_object_or_404(ApplicationUser, username=username).is_appuser == False:
            raise forms.ValidationError(
					self.error_messages['invalid_login'],
					code='user is not appuser',
					params={'username': self.username_field.verbose_name},
				)

        return self.cleaned_data


#=================================================================================================
# Application User
#=================================================================================================

#------------------------------------------------------------------------
Permission_Choices=[ ("Read","Read"),("Write","Write"),("Delete","Delete"), ("Admin","Admin"), ("No","No")]

class AppUser_Form(forms.ModelForm):

    permission=forms.ChoiceField(choices=Permission_Choices)

    class Meta:
        model=ApplicationUser
        fields=['first_name','last_name','email', 'is_active', 'username', 'name', 'initials','permission','is_appuser']

#------------------------------------------------------------------------
class AppUser_Filter(Filterbase_base):

    Search_all_fields = CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)

    username = CharFilter(lookup_expr='icontains',label='UQ Username')
    first_name = CharFilter(lookup_expr='icontains')
    last_name = CharFilter(lookup_expr='icontains')
    initials = CharFilter(lookup_expr='icontains')
    is_active = BooleanFilter(widget=forms.RadioSelect(choices=((True,'Yes'),(False,'No'))), label='Is Active')
    is_appuser = BooleanFilter(widget=forms.RadioSelect(choices=((True,'Yes'),(False,'No'))), label='Is AppUser')
    permission= ChoiceFilter(choices=Permission_Choices)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.form.initial['is_active'] = (True,'Yes')
        self.form.initial['is_appuser'] = (True,'Yes')

    class Meta:
        model=ApplicationUser
        fields=['username','first_name', 'last_name', 'initials','is_active','is_appuser','permission']

#=================================================================================================
# Application Dictionary
#=================================================================================================

class Dictionary_Form(forms.ModelForm):
    class Meta:
        model=Dictionary
        fields='__all__'

#------------------------------------------------------------------------
class Dictionary_Filter(Filterbase):
    dict_class = ChoiceFilter(choices=[])
    dict_value = CharFilter(lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        a=[tuple(d.values()) for d in Dictionary.objects.filter(astatus__gte=0).distinct().order_by('dict_class').values('dict_class')]
        choice_class=[(x[0], x[0]) for x in a]
        self.filters['dict_class'].extra["choices"] = choice_class
        self.filters['dict_class'].label='Class'
      
    class Meta:
        model=Dictionary
        fields=['dict_class']

#------------------------------------------------------------------------
## Document
class Document_Form(forms.ModelForm):
    doc_file = forms.FileField(label='Select a file', 
                                #   validators=[validate_file], 
                                  required=True)
    class Meta:
        model=Document
        fields='__all__'
   

#------------------------------------------------------------------------
## Image
#from .utils.form_wizard_tools import SelectFile_StepForm, MultipleFileField
#from .utils.files_upload import validate_file

# class Image_form(forms.ModelForm):
#     image_file = forms.ImageField(label='Select an image', 
#                                 #   validators=[validate_file], 
#                                   required=True)
    
#     class Meta:
#         model=Image
#         fields='__all__'


#===================================================================
# --Filterset Form--
#===================================================================

#------------------------------------------------------------------------
## Application User

    # @property
    # def qs(self):
    #     parent = super().qs
    #     print(parent)
    #     return parent.filter(is_appuser=True,is_active=True)
    
    # def filter_all_fields(self, queryset, name, value):
    #     if value:
    #         exclude_fields = ['password',]
    #         q_object = get_all_fields_q_object(self._meta.model, value,exclude_fields=exclude_fields)
    #         return queryset.filter(q_object)
    #     return queryset


#=================================================================================================
# Application log 
#=================================================================================================

class AppLog_Filter(Filterbase_base):
    Search_all_fields = CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)
    log_code = ChoiceFilter(choices=[])
    # start_date = DateFilter(field_name='log_time',lookup_expr=('gt'), widget=forms.DateInput(attrs={'type': 'date'}), label = 'Log date start from') 
    log_time = DateRangeFilter(label = 'Log Time')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['log_code'].extra["choices"] = self.Meta.model.get_field_choices(field_name='log_code')
      
    class Meta:
        model=ApplicationLog
        fields=list(model.HEADER_FIELDS.keys())


