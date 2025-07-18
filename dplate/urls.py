from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dplate.views import  (TestPlate_MapView, 
                    ) 
urlpatterns = [
    # Project 
    # path('project_card', Project_CardView.as_view(), name="project_card"),
    #path('project_list', Project_ListView.as_view(), name="project_list"),
    path('testplate/<str:pk>', TestPlate_MapView, name="testplate_map"),
    #path('createProject/', Project_CreateView, name="project_create"),
    #path('updateProject/<str:pk>', Project_UpdateView, name="project_update"),
    #path('deleteProject/<str:pk>', Project_RemoveView.as_view(), name="project_delete"),
    #path('project/report/<str:pk>', Project_ReportView, name="project_report"),

]
