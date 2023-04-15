"""
General base view class or function used by all applications
"""
import os
# from pathlib import Path
from datetime import datetime
from asgiref.sync import sync_to_async

from django import forms
from django.shortcuts import HttpResponse
from django.http import JsonResponse
# from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.generic import ListView
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
# from django.conf import settings


# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")

# ==========Super UserRequire Mixin===================================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = '/'

    def test_func(self):
        return self.request.user.has_permission('Admin')



# Model Validation utilities  =====================================================================================



    
# #-----------------------------------------------------------------------------------
# def instance_dict(instance, key_format=None):
#     "Returns a dictionary containing field names and values for the given instance"
# #-----------------------------------------------------------------------------------
#     from django.forms.models import model_to_dict
#     model_to_dict(instance, fields=[field.name for field in instance._meta.fields]) 



 
    


# ===========================================================================



#file path on server:

# Define full path name
# def get_filewithpath( file_name=None):
#     if settings.DEVELOPMENT=="Local":
#         file_path = f"static/images/{file_name}.svg"
   
#     else:
#         Base_dir = Path(__file__).resolve().parent.parent.parent
#         FILES_DIR=os.path.abspath(os.path.join(Base_dir, 'static/images'))
#         file_path=os.path.join(FILES_DIR, f"{file_name}.svg")
#     return file_path
# 
from django.urls import reverse_lazy
from django.views.generic.edit import DeleteView
class DeleteView_FKeyExist(DeleteView):
    model = None
    success_url = reverse_lazy("/")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["fkey_number"] = 20 # call fkey sum function
        return context

class ModelDeleteView(DeleteView):
    model = None
    success_url = reverse_lazy('/')
    template_name = None

    # def get_context_data(self, **kwargs):
    #     context = super().get_context_data(**kwargs)
    #     # context['can_delete'] = self.get_object().can_delete
    #     return context

    # def delete(self, request, *args, **kwargs):
    #     self.object = self.get_object()
    #     # if self.object.can_delete:
    #     if self.object:
    #         self.object.delete()
    #         data = {'message': 'Deleted successfully', 'deleted': True}
    #     else:
    #         data = {'message': 'Cannot delete the item', 'deleted': False}
    #     return JsonResponse(data)