
from django import forms
from .models import  ApplicationUser, Dictionaries


class ApplicationUser_form(forms.ModelForm):
    class Meta:
        model = ApplicationUser      
        fields= ['user_id', 'permissions', 'is_appuser']
        


# #=======================================Dictionary Form===========================================================
class Dictionary_form(forms.ModelForm):

    class Meta:
        model=Dictionaries
        fields="__all__"