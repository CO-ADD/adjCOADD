"""coadd_web URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

from app.views import index, userprofile, GroupListView, GroupCreateView, GroupUpdateView, GroupDeleteView, UserListView
from aa_chem import views

urlpatterns = [
    path('admin', admin.site.urls),
    path('', index, name="index"),
    path('compounds/', views.home, name="compounds"),
    path('usermanagement/', GroupListView.as_view(), name="usermanage"),
    path('usermanagement_create/', GroupCreateView.as_view(), name="createGroup"),
    path('usermanagement_update/<int:pk>', GroupUpdateView.as_view(), name="updateGroup"),
    path('usermanagement_delete/<int:pk>', GroupDeleteView.as_view(), name="deleteGroup"),
    path('userslist/', UserListView.as_view(), name="userslist"),

    path('accounts/', include('django.contrib.auth.urls')),
  
    path('importData/', views.importCSV, name="dataimport"),
    path('exportData/', views.exportCSV, name="dataexport"),
    path('accounts/', include('django.contrib.auth.urls')),
    path('user-profile/<int:id>/', userprofile, name='userprofile'),
    path('compounds/createTaxo', views.TaxoCreateView.as_view(), name="taxo_create"),
    path('compounds/createOrg', views.OrgCreateView.as_view(), name="org_create"),
    path('compounds/organism', views.OrgListView.as_view(), name="org_list"),
    path('compounds/taxo', views.TaxoListView.as_view(), name="taxo_list"),
    path('compounds/orgTable', views.OrgTableView.as_view(), name="org_table"),
    path('compounds/taxo/<int:pk>', views.TaxoUpdateView.as_view(), name="taxo_update"),
]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_URL)
    urlpatterns +=static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
