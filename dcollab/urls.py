from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dcollab.views import  (Organisation_ListView, Organisation_CreateView, Organisation_UpdateView,
                            CollabGroup_ListView, CollabGroup_CreateView,
                            CollabUser_ListView, CollabUser_CreateView,
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
    path('updateOrganisation/<str:pk>', Organisation_UpdateView.as_view(), name="organisation_update"),
    #path('deleteOrganisation/<str:pk>', Organisationt_RemoveView.as_view(), name="organisation_delete"),
    #path('organisation/report/<str:pk>', Organisation_ReportView, name="organisation_report"),

    #-- Collab Group 
    # path('collabgroup_card', CollabGroup_CardView.as_view(), name="collabgroup_card"),
    path('collabgroup_list', CollabGroup_ListView.as_view(), name="collabgroup_list"),
    #path('organisation/<str:pk>', CollabGroup_DetailView, name="collabgroup_detail"),
    path('createCollabGroup/', CollabGroup_CreateView, name="collabgroup_create"),
    #path('updateCollabGroup/<str:pk>', CollabGroup_UpdateView.as_view(), name="collabgroup_update"),
    #path('deleteCollabGroup/<str:pk>', CollabGroup_RemoveView.as_view(), name="collabgroup_delete"),
    #path('collabgroup/report/<str:pk>', CollabGroup_ReportView, name="collabgroup_report"),

    #-- Collab Group 
    # path('collabuser_card', CollabUser_CardView.as_view(), name="collabuser_card"),
    path('collabuser_list', CollabUser_ListView.as_view(), name="collabuser_list"),
    #path('organisation/<str:pk>', CollabUser_DetailView, name="collabuser_detail"),
    path('createCollabUser/', CollabUser_CreateView, name="collabuser_create"),
    #path('updateCollabUser/<str:pk>', CollabUser_UpdateView.as_view(), name="collabuser_update"),
    #path('deleteCollabUser/<str:pk>', CollabUser_RemoveView.as_view(), name="collabuser_delete"),
    #path('collabuser/report/<str:pk>', CollabUser_ReportView, name="collabuser_report"),
]
