
from django import forms
from .models import  ApplicationUser, Dictionaries


class ApplicationUser_form(forms.ModelForm):
    
    def save(self, *args, **kwargs):
        instance=super(ApplicationUser_form, self).save(commit=False)
        if instance.permission=='staff':
            instance.is_staff=True
            print('it working')
        instance.save() 
    class Meta:
        model = ApplicationUser      
        fields= ['user_id', 'permissions', 'is_appuser']
    
        


# #=======================================Dictionary Form===========================================================
class Dictionary_form(forms.ModelForm):

    class Meta:
        model=Dictionaries
        fields="__all__"