#import django_filters
from django import forms
from .models import  ApplicationUser, ApplicationLog, Dictionary, Document, Menu, Menu_Item
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate
from django.shortcuts import get_object_or_404
from django_filters import DateRangeFilter, CharFilter, BooleanFilter, ChoiceFilter, DateFromToRangeFilter, DateFilter

from applib.django.base.filters import BaseStatus_Filter, Base_Filter


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
PERMISSION_CHOICES=[ 
                    ("Read","Read"),
                    ("Write","Write"),
                    ("Delete","Delete"), 
                    ("Admin","Admin"), 
                    ("No","No")
                    ]

class AppUser_Form(forms.ModelForm):

    permission=forms.ChoiceField(choices=PERMISSION_CHOICES)

    class Meta:
        model=ApplicationUser
        fields=['first_name','last_name','email', 'is_active', 'username', 'name', 'initials','permission','is_appuser']

#------------------------------------------------------------------------
class AppUser_Filter(Base_Filter):

    Search_all_fields = CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)

    username = CharFilter(lookup_expr='icontains',label='UQ Username')
    first_name = CharFilter(lookup_expr='icontains')
    last_name = CharFilter(lookup_expr='icontains')
    initials = CharFilter(lookup_expr='icontains')
    is_active = BooleanFilter(widget=forms.RadioSelect(choices=((True,'Yes'),(False,'No'))), label='Is Active')
    is_appuser = BooleanFilter(widget=forms.RadioSelect(choices=((True,'Yes'),(False,'No'))), label='Is AppUser')
    permission= ChoiceFilter(choices=PERMISSION_CHOICES)

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


class Dictionary_UpdateForm(forms.ModelForm):
    dict_value = forms.CharField(disabled=True)

    class Meta:
        model = Dictionary
        fields = ['dict_value', 'dict_class', 'dict_desc', 'dict_sort', 'dict_app']

#------------------------------------------------------------------------
class Dictionary_Filter(BaseStatus_Filter):
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
   

#=================================================================================================
# Application log 
#=================================================================================================

class AppLog_Filter(Base_Filter):
    Search_all_fields = CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)
    log_code = ChoiceFilter(choices=[])
    # start_date = DateFilter(field_name='log_time',lookup_expr=('gt'), widget=forms.DateInput(attrs={'type': 'date'}), label = 'Log date start from') 
    log_time = DateRangeFilter(label = 'Log Time')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['log_code'].extra["choices"] = self.Meta.model.get_field_choices(field_name='log_code')
      
    class Meta:
        model=ApplicationLog
        fields=list(model.LIST_VIEW_FIELDS.keys())


#=================================================================================================
# Navigation Menu
#=================================================================================================
class Menu_Filter(BaseStatus_Filter):
    menu_name = CharFilter(field_name='menu_name', lookup_expr='icontains', label='Name')
    menu_label = CharFilter(field_name='menu_label', lookup_expr='icontains', label='Label')
    permission = CharFilter(field_name='permission', lookup_expr='icontains', label='Permission')

    class Meta:
        model = Menu
        fields = ['menu_name', 'menu_label', 'permission']


class Menu_Form(forms.ModelForm):
    class Meta:
        model = Menu
        exclude = ['menu_id']
        field_order = ['menu_name', 'menu_label', 'icon', 'order', 'permission']


class Menu_UpdateForm(Menu_Form):
    menu_id = forms.CharField(disabled=True)

    class Meta:
        model = Menu
        exclude = []
        field_order = ['menu_id', 'menu_name', 'menu_label', 'icon', 'order', 'permission']


# -----------------------------------------------------------------
class Menu_Item_Filter(BaseStatus_Filter):
    label = CharFilter(field_name='label', lookup_expr='icontains', label='Label')
    url_name = CharFilter(field_name='url_name', lookup_expr='icontains', label='URL Name')
    permission = CharFilter(field_name='permission', lookup_expr='icontains', label='Permission')
    menu_id = ChoiceFilter(field_name='menu_id', choices=[], label='Menu')
    parent = ChoiceFilter(field_name='parent', choices=[], label='Parent')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['menu_id'].extra['choices'] = [
            (m.menu_id, f"[{m.menu_id}] {m.menu_name}") for m in Menu.objects.all().order_by('order', 'menu_name')
        ]
        self.filters['parent'].extra['choices'] = [
            (i.item_id, f"[{i.item_id}] {i.label}") for i in Menu_Item.objects.select_related('menu_id').all().order_by('menu_id', 'order', 'label')
        ]

    class Meta:
        model = Menu_Item
        fields = ['label', 'url_name', 'permission', 'menu_id', 'parent']


class Menu_Item_Form(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        super(Menu_Item_Form, self).__init__(*args, **kwargs)
        for field_name in self.fields:
            self.fields[field_name].label = self.Meta.model._meta.get_field(field_name).verbose_name

    class Meta:
        model = Menu_Item
        exclude = ['item_id']
        field_order = ['menu_id', 'parent', 'label', 'url_name', 'url_path', 'icon', 'order',
                       'separator', 'disabled', 'modal', 'permission']


class Menu_Item_UpdateForm(Menu_Item_Form):
    item_id = forms.CharField(disabled=True)

    class Meta:
        model = Menu_Item
        exclude = []
        field_order = ['item_id', 'menu_id', 'parent', 'label', 'url_name', 'url_path', 'icon', 'order',
                       'separator', 'disabled', 'modal', 'permission']

    