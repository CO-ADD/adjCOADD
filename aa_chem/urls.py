from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include

from aa_chem.views import  (home, exportCSV, import_excel_taxo,import_excel_dict, createOrgnisms, detailOrganism,
                        updateOrganism, TaxoListView,TaxoCreate,TaxoUpdate,deleteTaxonomy, OrganismListView, OrganismCardView, deleteOrganism)
from aa_chem.utils import searchbar_01


urlpatterns = [

    path('', home, name="compounds"),
    path('taxo', TaxoListView.as_view(), name="taxo_list"),
    path('taxoListview', TaxoListView.as_view(), name="taxo_list_view"),
    path('organism_list', OrganismListView.as_view(), name="org_list"),
    path('organism_card', OrganismCardView.as_view(), name="org_card"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
    

    path('searchbar_01/', searchbar_01, name="searchbar_01"),
    path('createOrg/', createOrgnisms, name="org_create"),
    path('createTaxo/', TaxoCreate, name="taxo_create"),

    path('taxo/<str:pk>', TaxoUpdate, name="taxo_update"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),

    path('deleteOrg/<str:pk>', deleteOrganism, name="organism_delete"),
    path('deleteTaxo/<str:pk>', deleteTaxonomy, name="taxonomy_delete"),
    #=======================Data Export/Import===================================================================
    path('exportData/', exportCSV, name="dataexport"),
    path('import_Taxonomy/', import_excel_taxo, name="importTaxo"),
    path('import_dictionary/', import_excel_dict, name="importDict"),
]