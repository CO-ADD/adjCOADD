from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

from app.views import (index, userprofile, AppUserListView, AppUserCreateView, AppUserUpdateView, 
    AppUserDeleteView, AppUserListView, DictionariesView, createDictionary, deleteDictionary)


urlpatterns = [

    path('user_list/', AppUserListView.as_view(), name="userslist"),
    path('user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('user_update/<int:pk>', AppUserUpdateView.as_view(), name="updateAppUser"),
    path('user_delete/<int:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('dict/', DictionariesView.as_view(), name='dict_view' ),
    path('dict/create', createDictionary, name='dict_create' ),
    path('dictionary_update/<str:pk>', createDictionary, name='updateDictionary' ),

    path('dictionary_delete/<str:pk>', deleteDictionary, name="deleteDictionary"),
]