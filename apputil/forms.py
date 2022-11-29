
from django import forms
from .models import  ApplicationUser, Dictionary
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth import authenticate


class ApplicationUser_form(forms.ModelForm):
    pass
    
#     def save(self, *args, **kwargs):
#         instance=super(ApplicationUser_form, self).save(commit=False)
#         if instance.permission=='staff':
#             instance.is_staff=True
#             print('it working')
#         instance.save() 
#     class Meta:
#         model = ApplicationUser      
#         fields= ['user_id', 'permissions', 'is_appuser']
    
        


# #=======================================Dictionary Form===========================================================
class Dictionary_form(forms.ModelForm):
    class Meta:
        model=Dictionary
        fields='__all__'


class Login_form(AuthenticationForm):
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

