from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import (index, userprofile, AppUserListView, AppUserCreateView, updateApplicationuser, 
    AppUserDeleteView, AppUserListView, DictionaryView, createDictionary,updateDictionary, deleteDictionary,exportCSV, testsite) 



urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUserListView.as_view(), name="userslist"),
    path('user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', updateApplicationuser, name="updateAppUser"),
    path('user_delete/<str:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('dict/', DictionaryView.as_view(), name='dict_view' ),
    path('dict/create', createDictionary, name='dict_create' ),
    path('dict_update', updateDictionary, name='dict_update' ),
    path('dict_delete', deleteDictionary, name='dict_delete' ),
    #Data Export/Import 
    path('exportData/', exportCSV, name="dataexport"),
    path('testsite/', testsite, name="testsite"),

   
]