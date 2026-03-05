from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dsample.views import  (Project_ListView, Project_CreateView, Project_DetailView, Project_UpdateView, Project_ReportView,
                            Project_StockPrepView,
                            # Project_RemoveView,
                            # ScreenRun_CreateView, ScreenRun_UpdateView, 
                            # Assay_ListView,Assay_DetailView,Assay_CreateView, Assay_UpdateView, Assay_RemoveView,
                    )
from dsample.process_views import (Load_StockPrep_ProcessView, Load_Project_ProcessView, Add_ProjectInfo_ProcessView)
 
urlpatterns = [
    # Project 
    # path('project_card', Project_CardView.as_view(), name="project_card"),
    path('project/list', Project_ListView.as_view(), name="project_list"),
    path('project/<str:pk>', Project_DetailView, name="project_detail"),
    path('Project/create', Project_CreateView, name="project_create"),
    path('project/update/<str:pk>', Project_UpdateView, name="project_update"),
    #path('deleteProject/<str:pk>', Project_RemoveView.as_view(), name="project_delete"),
    path('project/report/<str:pk>', Project_ReportView, name="project_report"),

    path('Project/load_submission', Load_Project_ProcessView.as_view(), name='load_submission'),
    path('project/add_projectinfo/<str:pk>', Add_ProjectInfo_ProcessView.as_view(), name="add_projectinfo"),
    path('project/generate_stockprep/<str:pk>', Project_StockPrepView, name="generate_stockprep"),
    path('project/load_stockprep/<str:pk>', Load_StockPrep_ProcessView.as_view(), name='load_stockprep'),
]
