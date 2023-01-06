from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (DrugListView, DrugCardView, detailDrug, createDrug, updateDrug, deleteDrug)



urlpatterns = [
    path('drug_card', DrugCardView.as_view(), name="drug_card"),
    path('drug_list', DrugListView.as_view(), name="drug_list"),
    path('drug/<str:pk>', detailDrug, name="drug_detail"),
    path('createDrug/', createDrug, name="drug_create"),
    path('updateDrug/<str:pk>', updateDrug, name="drug_update"),
    path('deleteDrug/<str:pk>', deleteDrug, name="drug_delete"),
]