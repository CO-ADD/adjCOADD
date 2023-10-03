from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from apputil.views import (index, userprofile, 
                           AppUser_ListView, AppUser_CreateView, AppUser_UpdateView,  AppUser_DetailView, AppUser_DeleteView, 
                           AppLog_ListView, 
                           Dictionary_ListView, Dictionary_CreateView,updateDictionary, deleteDictionary,
                            DataExportView, Importhandler_apputils, CreatedocumentView, DocDeleteView)

from apputil.utils.flex_pivottable import flex_pivottable




urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUser_ListView.as_view(), name="userslist"),
    path('user_create/', AppUser_CreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', AppUser_UpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<str:pk>', AppUser_DeleteView.as_view(), name="deleteAppUser"),
    path('user_profile/<str:pk>', AppUser_DetailView.as_view(), name='userprofile' ),
    
    path('log_list/', AppLog_ListView.as_view(), name='loglist' ),
    
    path('dict/', Dictionary_ListView.as_view(), name='dict_view' ),
    path('dict_create/', Dictionary_CreateView.as_view(), name='dict_create' ),
    path('dict_update/', updateDictionary, name='dict_update' ),
    path('dict_delete/', deleteDictionary, name='dict_delete' ),
 
    # path('img/<str:pk>', CreateimageView.as_view(), name="addimg"),
    path('doc/<str:pk>', CreatedocumentView.as_view(), name="adddoc"),
    # path('img-delete/<str:pk>', ImageDeleteView.as_view(), name='org_img_delete'),
    path('doc-delete/<str:pk>', DocDeleteView.as_view(), name='org_doc_delete'),
    # path('data-visual/<str:process_name>', Data_visualView.as_view(), name="data-visual"),
    
    path('exportData/', DataExportView.as_view(), name="dataexport"),
    path('import-excel/<str:process_name>', Importhandler_apputils.as_view(), name="excel-import"),
    path('pivotedtableview/<str:app_model>',flex_pivottable, name="pivoted-table"),
]