
from django import forms
from .models import Groupfilter, User


class GroupCreate(forms.ModelForm):
    class Meta:
        model = Groupfilter
       
        fields='__all__'
        widgets = {
            'text': forms.TextInput(attrs={
                'id': 'post-text', 
                'required': True, 
                'placeholder': 'Say something...'
            }),
        }