from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from apputil.views import (index, userprofile, 
                           AppUser_ListView, AppUser_CreateView, AppUser_UpdateView,  AppUser_DetailView, AppUser_RemoveView, 
                           AppLog_ListView, 
                           Dictionary_ListView, Dictionary_CreateView, Dictionary_UpdateView, updateDictionary, removeDictionary,
                           DataExportView, Importhandler_apputils, CreatedocumentView, DocDeleteView,
                           SystemInfoView,
                           Menu_ListView, Menu_CreateView, Menu_UpdateView,
                           Menu_Item_ListView, Menu_Item_CreateView, Menu_Item_UpdateView)

from apputil.utils.flex_pivottable import flex_pivottable
from apputil.test_app import Test_View
from applib.process.process_forms import  Start_Progress

from apputil.task_views import Task_UpdateABaseView, TaskResult_ListView, Task_RemoveView

urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUser_ListView.as_view(), name="userslist"),
    path('user_create/', AppUser_CreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', AppUser_UpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<str:pk>', AppUser_RemoveView.as_view(), name="deleteAppUser"),
    path('user_profile/<str:pk>', AppUser_DetailView.as_view(), name='userprofile' ),
    
    path('log_list/', AppLog_ListView.as_view(), name='loglist' ),
    
    path('dict/', Dictionary_ListView.as_view(), name='dict_view' ),
    path('dict_create/', Dictionary_CreateView.as_view(), name='dict_create' ),
    path('dict_update/<str:pk>', Dictionary_UpdateView.as_view(), name='dict_update' ),
    path('dict_delete/', removeDictionary, name='dict_delete' ),

    # Navigation Menu
    path('menu_list/', Menu_ListView.as_view(), name='menu_list'),
    path('menu_create/', Menu_CreateView, name='menu_create'),
    path('menu_update/<str:pk>', Menu_UpdateView.as_view(), name='menu_update'),
    path('menuitem_list/', Menu_Item_ListView.as_view(), name='menuitem_list'),
    path('menuitem_create/', Menu_Item_CreateView, name='menuitem_create'),
    path('menuitem_update/<str:pk>', Menu_Item_UpdateView.as_view(), name='menuitem_update'),
 
    # path('img/<str:pk>', CreateimageView.as_view(), name="addimg"),
    path('doc/<str:pk>', CreatedocumentView.as_view(), name="adddoc"),
    # path('img-delete/<str:pk>', ImageDeleteView.as_view(), name='org_img_delete'),
    path('doc-delete/<str:pk>', DocDeleteView.as_view(), name='org_doc_delete'),
    # path('data-visual/<str:process_name>', Data_visualView.as_view(), name="data-visual"),

    path('progress/', Start_Progress, name='start_progress'),
    path('info/', SystemInfoView, name='system_info'),

    # Background tasks (django.tasks)
    path('tasks/update-abase/', Task_UpdateABaseView, name='task_update_abase'),
    path('tasks/current/', TaskResult_ListView.as_view(), name='task_current'),
    path('tasks/remove/<str:pk>', Task_RemoveView, name='task_remove'),
    
    path('exportData/', DataExportView.as_view(), name="dataexport"),
    path('import-excel/<str:process_name>', Importhandler_apputils.as_view(), name="excel-import"),
    path('pivotedtableview/<str:app_model>',flex_pivottable, name="pivoted-table"),
    
    path('test',Test_View,name='test-view'),
]