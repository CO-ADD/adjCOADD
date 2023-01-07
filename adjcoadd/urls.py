"""adjCOADD URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
"""
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

import dorganism.urls
import apputil.urls

from apputil.views import login_user, logout_user, permission_not_granted

urlpatterns = [
    path('admin/', admin.site.urls), 
    #========================appusers URL==========================================================================
    path('', login_user, name='login' ),
    path('accounts/login/', login_user, name='login'),
    path('accounts/logout/', logout_user, name='logout' ), 
    path('permission_not_granted/', permission_not_granted, name='permission_not_granted'),
    path('', include('apputil.urls')),
    #========================OrgDB model views URL====View, Create, Updata, Delete=================================
    path('', include('dorganism.urls')),
    # path('', include('ddrug.urls')),
]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_URL)
    urlpatterns +=static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
