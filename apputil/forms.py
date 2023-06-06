import django_filters
from django import forms
from .models import  ApplicationUser, Dictionary, Image, Document
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate
from django.shortcuts import get_object_or_404

from .utils.filters_base import Filterbase


# --Login Form--
class Login_form(AuthenticationForm):
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


# --Model Form--
## Application User

Permission_Choices=[ ("Read","Read"),("Write","Write"),("Delete","Delete"), ("Admin","Admin"), ("No","No")]

class ApplicationUser_form(forms.ModelForm):
    permission=forms.ChoiceField(choices=Permission_Choices)
    class Meta:
        model=ApplicationUser
        fields=['first_name','last_name','email', 'is_active', 'username', 'name', 'initials','permission','is_appuser']


## Dictionary
class Dictionary_form(forms.ModelForm):
    class Meta:
        model=Dictionary
        fields='__all__'

## Image
from .utils.form_wizard_tools import SelectFile_StepForm, MultipleFileField
from .utils.files_upload import validate_file

class Image_form(forms.ModelForm):
    image_file = forms.ImageField(label='Select an image', 
                                  validators=[validate_file], 
                                  required=False)
    

    class Meta:
        model=Image
        fields='__all__'

    

# --Filterset Form--
## Application User

class AppUserfilter(django_filters.FilterSet):
    Search_all_fields = django_filters.CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)

    username = django_filters.CharFilter(lookup_expr='icontains')
    first_name = django_filters.CharFilter(lookup_expr='icontains')
    last_name = django_filters.CharFilter(lookup_expr='icontains')
    permission=django_filters.ChoiceFilter(choices=Permission_Choices)
         
    class Meta:
        model=ApplicationUser
        fields=['username','first_name', 'last_name', 'permission']

    @property
    def qs(self):
        parent = super().qs
        return parent.filter(is_appuser=True)
    
    def filter_all_fields(self, queryset, name, value):
        if value:
            exclude_fields = ['password',]
            q_object = get_all_fields_q_object(self._meta.model, value,exclude_fields=exclude_fields)
            return queryset.filter(q_object)
        return queryset

## Dictionary
class Dictionaryfilter(Filterbase):
    dict_class = django_filters.ChoiceFilter(choices=[])
    dict_value = django_filters.CharFilter(lookup_expr='icontains')
    #   dict_desc = django_filters.CharFilter(lookup_expr='icontains')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        a=[tuple(d.values()) for d in Dictionary.objects.filter(astatus__gte=0).order_by().values('dict_class').distinct()]
        choice_class=[(x[0], x[0]) for x in a]
        self.filters['dict_class'].extra["choices"] = choice_class
        self.filters['dict_class'].label='Class'
      
    class Meta:
        model=Dictionary
        fields=['dict_class']


