
from django import forms
from .models import  ApplicationUser, Dictionaries


class GroupCreate(forms.ModelForm):
    class Meta:
        model = ApplicationUser
       
        fields='__all__'
        widgets = {
            'text': forms.TextInput(attrs={
                'id': 'post-text', 
                'required': True, 
                'placeholder': 'Say something...'
            }),
        }


# #=======================================Dictionary Form===========================================================
class Dictionary_form(forms.ModelForm):

    class Meta:
        model=Dictionaries
        fields="__all__"