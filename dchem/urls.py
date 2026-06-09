from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dchem.views import Compound_ListAPI


urlpatterns = [
    # Chem Structure
    path('api/chemstructure', Compound_ListAPI.as_view({'get': 'list'}), name="chemstructure_api_list"),
    #path('api/chemstructure/<str:pk>', Compound_DetailAPI.as_view({'get': 'retrieve',"patch": "partial_update","post": "update"}), name="chemstructure_api_detail"),
]
