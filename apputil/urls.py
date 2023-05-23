from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from apputil.views import (index, userprofile, AppUserListView, AppUserCreateView, ApplicationUserUpdateView, 
    AppUserDeleteView, AppUserListView, DictionaryView, DictionaryCreateView,updateDictionary, deleteDictionary,
    exportCSV, Importhandler_apputils) 



urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUserListView.as_view(), name="userslist"),
    path('user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', ApplicationUserUpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<str:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('dict/', DictionaryView.as_view(), name='dict_view' ),
    path('dict_create/', DictionaryCreateView.as_view(), name='dict_create' ),
    path('dict_update/', updateDictionary, name='dict_update' ),
    path('dict_delete/', deleteDictionary, name='dict_delete' ),
    path('exportData/', exportCSV, name="dataexport"),
    path('import-excel/<str:process_name>', Importhandler_apputils.as_view(), name="excel-import"),
    

   
]