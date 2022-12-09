from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import (index, userprofile, AppUserListView, AppUserCreateView, AppUserUpdateView, 
    AppUserDeleteView, AppUserListView, DictionaryView, createDictionary, Importhandler)


urlpatterns = [
    path('index/', index, name="index"),
    path('user_list/', AppUserListView.as_view(), name="userslist"),
    path('user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('user_update/<str:pk>', AppUserUpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<str:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('dict/', DictionaryView.as_view(), name='dict_view' ),
    path('dict/create', createDictionary, name='dict_create' ),
    path("import/", Importhandler.as_view(), name="import"),
    path("tasks/save", Importhandler.save_task, name="save_task"),
    path("tasks/<task_id>/", Importhandler.get_status, name="status_task"),
    path("tasks/", Importhandler.run_task, name="run_task"),
    path("tasks/cancel", Importhandler.delete_task, name="cancel_task"),
    # path("tasks/save", Importhandler.save_task, name="save_task"),

    # path('dictionary_update/<str:pk>', createDictionary, name='updateDictionary' ),
    # re_path(r'dictionary_update/(?P<pk>\w+)/', createDictionary, name='updateDictionary'),
    # re_path(r'dictionary_update/(?P<pk>\w+)/', deleteDictionary, name="deleteDictionary"),
    # path('dictionary_delete/<str:pk>', deleteDictionary, name="deleteDictionary"),
]