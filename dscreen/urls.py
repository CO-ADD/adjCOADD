from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dscreen.views import  (ScreenRun_ListView,ScreenRun_DetailView,ScreenRun_DeleteView,
                            # ScreenRun_CreateView, ScreenRun_UpdateView, 
                            # Assay_ListView,Assay_DetailView,Assay_CreateView, Assay_UpdateView, Assay_DeleteView,
                    ) 

urlpatterns = [
    # ScreenRun 
    # path('screenrun_card', ScreenRun_CardView.as_view(), name="screenrun_card"),
    path('screenrun_list', ScreenRun_ListView.as_view(), name="screenrun_list"),
    # path('screenrun/<str:pk>', ScreenRun_DetailView, name="screenrun_detail"),
    # path('createScreenrun/', ScreenRun_CreateView, name="screenrun_create"),
    # path('updateScreenrun/<str:pk>', ScreenRun_UpdateView, name="screenrun_update"),
    path('deleteScreenrun/<str:pk>', ScreenRun_DeleteView.as_view(), name="screenrun_delete"),

]
