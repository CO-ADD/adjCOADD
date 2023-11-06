from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (*) 
from .utils.utils import search_organism, search_organism_id


urlpatterns = [
    # modelpath
    # path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
  
]