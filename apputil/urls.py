from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from apputil.views import (index, userprofile, AppUserListView, AppUserCreateView, ApplicationUserUpdateView,  AppUserDetailView,
    AppUserDeleteView, AppUserListView, AppLogView, DictionaryView, DictionaryCreateView,updateDictionary, deleteDictionary,
    DataExportView, Importhandler_apputils, CreatedocumentView, DocDeleteView)

from .utils.flex_pivottable import flex_pivottable




urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUserListView.as_view(), name="userslist"),
    path('user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', ApplicationUserUpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<str:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('user_profile/<str:pk>', AppUserDetailView.as_view(), name='userprofile' ),
    path('log_list/', AppLogView.as_view(), name='loglist' ),
    path('dict/', DictionaryView.as_view(), name='dict_view' ),
    path('dict_create/', DictionaryCreateView.as_view(), name='dict_create' ),
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