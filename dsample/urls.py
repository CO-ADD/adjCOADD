from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dsample.views import  (Project_ListView, Project_CreateView, Project_DetailView, Project_UpdateView, Project_ReportView,
                            Project_CmpdPrepView,
                            # Project_RemoveView,
                            # ScreenRun_CreateView, ScreenRun_UpdateView, 
                            # Assay_ListView,Assay_DetailView,Assay_CreateView, Assay_UpdateView, Assay_RemoveView,
                    ) 
urlpatterns = [
    # Project 
    # path('project_card', Project_CardView.as_view(), name="project_card"),
    path('project/list', Project_ListView.as_view(), name="project_list"),
    path('project/<str:pk>', Project_DetailView, name="project_detail"),
    path('Project/create', Project_CreateView, name="project_create"),
    path('Project/<str:pk>/edit', Project_UpdateView, name="project_update"),
    #path('deleteProject/<str:pk>', Project_RemoveView.as_view(), name="project_delete"),
    path('project/<str:pk>/report', Project_ReportView, name="project_report"),

    path('project/<str:pk>/generate_cmpdprep', Project_CmpdPrepView, name="generate_cmpdprep"),
]
