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

from app.views import index, userprofile, AppUserListView, AppUserCreateView, AppUserUpdateView, AppUserDeleteView, AppUserListView, DictionariesView,DictCreate
from aa_chem.views import  (home, exportCSV, import_excel_taxo,import_excel_dict, searchbar_01, createOrgnisms, organismDetail,
                        updateOrganism, TaxoListView,TaxoCreate,TaxoUpdate,deleteTaxonomy, OrgListView, OrgCardView, deleteOrganism)

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', index, name="index"),

    #========================appusers URL==========================================================================
    path('accounts/', include('django.contrib.auth.urls')),  
    path('app/user_list/', AppUserListView.as_view(), name="userslist"),
    path('app/user_create/', AppUserCreateView.as_view(), name="createAppUser"),
    path('app/user_update/<int:pk>', AppUserUpdateView.as_view(), name="updateAppUser"),
    path('app/user_delete/<int:pk>', AppUserDeleteView.as_view(), name="deleteAppUser"),
    path('app/dict/', DictionariesView.as_view(), name='dict_view' ),
    path('app/dict/create', DictCreate, name='dict_create' ),

    #========================OrgDB model views URL====View, Create, Updata, Delete=================================
    path('aa_chem/', home, name="compounds"),
    path('aa_chem/taxo', TaxoListView.as_view(), name="taxo_list"),
    path('aa_chem/taxoListview', TaxoListView.as_view(), name="taxo_list_view"),
    path('aa_chem/organism_list', OrgListView.as_view(), name="org_list"),
    path('aa_chem/organism_card', OrgCardView.as_view(), name="org_card"),
    path('aa_chem/organism/<str:pk>', organismDetail, name="org_detail"),
    

    path('aa_chem/searchbar_01/', searchbar_01, name="searchbar_01"),
    path('aa_chem/createOrg/', createOrgnisms, name="org_create"),
    path('aa_chem/createTaxo/', TaxoCreate, name="taxo_create"),

    path('aa_chem/taxo/<str:pk>', TaxoUpdate, name="taxo_update"),
    path('aa_chem/updateOrg/<str:pk>', updateOrganism, name="organism_update"),

    path('aa_chem/deleteOrg/<str:pk>', deleteOrganism, name="organism_delete"),
    path('aa_chem/deleteTaxo/<str:pk>', deleteTaxonomy, name="taxonomy_delete"),
    #========================Data Export/Import===================================================================
    path('aa_chem/exportData/', exportCSV, name="dataexport"),
    path('aa_chem/import_Taxonomy/', import_excel_taxo, name="importTaxo"),
    path('aa_chem/import_dictionary/', import_excel_dict, name="importDict"),
]

if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_URL)
    urlpatterns +=static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
