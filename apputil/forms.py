from django import forms
from .models import  ApplicationUser, Dictionary
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate


class ApplicationUser_form(forms.ModelForm):
    pass
    

# #=======================================Dictionary Form===========================================================
class Dictionary_form(forms.ModelForm):
    class Meta:
        model=Dictionary
        fields='__all__'


class Login_form(AuthenticationForm):
    username=forms.CharField(widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'UQ Username', }), label='')
    password= forms.CharField(widget=forms.PasswordInput(
    attrs={'class':'form-control','type':'password', 'name': 'password','placeholder':'Password'}),
    label='')
    def clean(self):
        username = self.cleaned_data.get('username')
        password = self.cleaned_data.get('password')
        print(username)
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
       

        return self.cleaned_data

# ==============Validators=======================================================
# import magic


from django.utils.deconstruct import deconstructible
from django.template.defaultfilters import filesizeformat
from django.core.exceptions import ValidationError


@deconstructible
class FileValidator(object):
    error_messages = {
     'max_size': ("Ensure this file size is not greater than %(max_size)s."
                  " Your file size is %(size)s."),
     'min_size': ("Ensure this file size is not less than %(min_size)s. "
                  "Your file size is %(size)s."),
     'content_type': "Files of type %(content_type)s are not supported.",
    }

    def __init__(self, max_size=None, min_size=None, content_types=()):
        self.max_size = max_size
        self.min_size = min_size
        self.content_types = content_types

    def __call__(self, data):
        if self.max_size is not None and data.size > self.max_size:
            params = {
                'max_size': filesizeformat(self.max_size), 
                'size': filesizeformat(data.size),
            }
            raise ValidationError(self.error_messages['max_size'],
                                   'max_size', params)

        if self.min_size is not None and data.size < self.min_size:
            params = {
                'min_size': filesizeformat(self.min_size),
                'size': filesizeformat(data.size)
            }
            raise ValidationError(self.error_messages['min_size'], 
                                   'min_size', params)

        # if self.content_types:
        #     content_type = magic.from_buffer(data.read(), mime=True)
        #     data.seek(0)

        #     if content_type not in self.content_types:
        #         params = { 'content_type': content_type }
        #         raise ValidationError(self.error_messages['content_type'],
        #                            'content_type', params)

    def __eq__(self, other):
        return (
            isinstance(other, FileValidator) and
            self.max_size == other.max_size and
            self.min_size == other.min_size and
            self.content_types == other.content_types
        )