from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (DrugListView, DrugCardView,detailDrug, createDrug, updateDrug, 
    detailVitekcard,  VitekcardListView, Importhandler_VITEK, smartsQuery)#VitekcardListView,



urlpatterns = [
    path('drug_card', DrugCardView.as_view(), name="drug_card"),
    path('drug_list', DrugListView.as_view(), name="drug_list"),
    path('drug/<str:pk>', detailDrug, name="drug_detail"),
    path('drug_detail_structure/<str:pk>', smartsQuery, name="smartsquery"),
    path('createDrug/', createDrug, name="drug_create"),
    path('updateDrug/<str:pk>', updateDrug, name="drug_update"),
    path('vitekcard_list', VitekcardListView.as_view(), name="vitekcard_list"),
    path('vitekcard_detail/<str:pk>', detailVitekcard, name="vitekcard_detail"),
    path("import-VITEK/", Importhandler_VITEK.as_view(), name="import-VITEK"),
]