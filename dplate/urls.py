from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dplate.views import (TestPlate_ListView, TestPlate_DetailView, TestPlate_UpdateView, TestPlate_MapView,
                          MasterPlate_ListView, MasterPlate_DetailView, MasterPlate_UpdateView,
                          Labware_ListView, Labware_UpdateView)

urlpatterns = [
    # TestPlate
    path('testplate_list', TestPlate_ListView.as_view(), name="testplate_list"),
    path('testplate/update/<str:pk>', TestPlate_UpdateView, name="testplate_update"),
    path('testplate/map/<str:pk>', TestPlate_MapView, name="testplate_map"),
    path('testplate/<str:pk>', TestPlate_DetailView, name="testplate_detail"),

    # MasterPlate
    path('masterplate_list', MasterPlate_ListView.as_view(), name="masterplate_list"),
    path('masterplate/update/<str:pk>', MasterPlate_UpdateView, name="masterplate_update"),
    path('masterplate/<str:pk>', MasterPlate_DetailView, name="masterplate_detail"),

    # Labware
    path('labware_list', Labware_ListView.as_view(), name="labware_list"),
    path('labware/update/<str:pk>', Labware_UpdateView, name="labware_update"),
]
