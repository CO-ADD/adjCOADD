from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dcollab.views import  (Organisation_ListView, Organisation_CreateView,
                            # Project_CreateView, Project_DetailView, Project_UpdateView, Project_ReportView,
                            # Project_RemoveView,
                            # ScreenRun_CreateView, ScreenRun_UpdateView, 
                            # Assay_ListView,Assay_DetailView,Assay_CreateView, Assay_UpdateView, Assay_RemoveView,
                    ) 
urlpatterns = [
    #-- Organisations 
    # path('organisation_card', Organisation_CardView.as_view(), name="organisation_card"),
    path('organisation_list', Organisation_ListView.as_view(), name="organisation_list"),
    #path('organisation/<str:pk>', Organisation_DetailView, name="organisation_detail"),
    path('createOrganisation/', Organisation_CreateView, name="organisation_create"),
    #path('updateOrganisation/<str:pk>', Organisation_UpdateView, name="organisation_update"),
    #path('deleteOrganisation/<str:pk>', Organisationt_RemoveView.as_view(), name="organisation_delete"),
    #path('organisation/report/<str:pk>', Organisation_ReportView, name="organisation_report"),

    #-- Organisations 

    #-- Organisations 
]
