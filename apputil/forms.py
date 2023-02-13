from django import forms
from .models import  ApplicationUser, Dictionary
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate

Permission_Choices=[ ("Read","Read"),("Write","Write"),("Delete","Delete"), ("Admin","Admin"), ("No","No")]

class ApplicationUser_form(forms.ModelForm):
    permission=forms.ChoiceField(choices=Permission_Choices)
    class Meta:
        model=ApplicationUser
        fields=['first_name','last_name','email', 'is_active', 'username', 'name', 'initials','permission','is_appuser']

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

