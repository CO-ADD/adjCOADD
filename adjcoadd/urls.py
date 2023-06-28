"""adjCOADD URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
"""
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include
from django.conf.urls import (
handler400, handler403, handler404, handler500
)
import dorganism.urls
import apputil.urls

from apputil.views import login_user, logout_user, permission_not_granted


handler404 = "apputil.views.custom_page_not_found_view"

urlpatterns = [
    path('admin/', admin.site.urls), 
    #========================appusers URL==========================================================================
    path('', login_user, name='login' ),
    path('accounts/login/', login_user, name='login'),
    path('accounts/logout/', logout_user, name='logout' ), 
    path('permission_not_granted/', permission_not_granted, name='permission_not_granted'),
    path('', include('apputil.urls')),
    #========================OrgDB model views URL====View, Create, Updata, Delete=================================
    path('dorganism/', include('dorganism.urls')),
    path('ddrug/', include('ddrug.urls')),
    path('dgene/', include('dgene.urls')),
]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns +=static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
