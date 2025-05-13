from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dscreen.views import  (ScreenRun_ListView,ScreenRun_CreateView,ScreenRun_DetailView,ScreenRun_UpdateView,ScreenRun_DeleteView,
                            # ScreenRun_CreateView, ScreenRun_UpdateView, 
                            # Assay_ListView,Assay_DetailView,Assay_CreateView, Assay_UpdateView, Assay_DeleteView,
                    ) 

urlpatterns = [
    # ScreenRun 
    # path('screenrun_card', ScreenRun_CardView.as_view(), name="screenrun_card"),
    path('screenrun_list', ScreenRun_ListView.as_view(), name="screenrun_list"),
    path('screenrun/<str:pk>', ScreenRun_DetailView, name="screenrun_detail"),
    path('createScreenrun/', ScreenRun_CreateView.as_view(), name="screenrun_create"),
    path('updateScreenrun/<str:pk>', ScreenRun_UpdateView.as_view(), name="screenrun_update"),
    path('deleteScreenrun/<str:pk>', ScreenRun_DeleteView.as_view(), name="screenrun_delete"),

]


    # path('drug_card', DrugCardView.as_view(), name="drug_card"),
    # path('drug_list', DrugListView.as_view(), name="drug_list"),
    # path('drug/<str:pk>', detailDrug, name="drug_detail"),
    # path('drug_detail_structure/<str:pk>', smartsQuery, name="smartsquery"),
    # path('createDrug/', DrugCreateView.as_view(), name="drug_create"),
    # path('updateDrug/<str:pk>', DrugUpdateView.as_view(), name="drug_update"),
