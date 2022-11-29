from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (exportCSV, import_excel_taxo,import_excel_dict, import_excel_organism, createOrganism, detailOrganism,
                        updateOrganism, TaxonomyCardView,TaxonomyListView,detailTaxonomy, createTaxonomy,updateTaxonomy, deleteTaxonomy, 
                       OrganismListView,OrganismCardView, deleteOrganism) # detailChangeOrganism, 
from .utils import searchbar_01


urlpatterns = [
    path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
    path('taxonomy_list', TaxonomyListView.as_view(), name="taxo_list"),
    path('taxonomy/<str:pk>', detailTaxonomy, name="taxo_detail"),
    path('organism_card', OrganismCardView.as_view(), name="org_card"),
    path('organism_list', OrganismListView.as_view(), name="org_list"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
   

    path('searchbar_01/', searchbar_01, name="searchbar_01"),
    path('createOrg/', createOrganism, name="org_create"),
    path('createTaxo/', createTaxonomy, name="taxo_create"),

    path('updatetaxo/<str:pk>', updateTaxonomy, name="taxonomy_update"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),
    # path('updateOrgdetail/', detailChangeOrganism, name="organism_updatedetail"),


    path('deleteOrg/<str:pk>', deleteOrganism, name="organism_delete"),
    re_path(r'^deleteTaxo/(?P<pk>\w+\s*\w+)/', deleteTaxonomy, name="taxonomy_delete"),
  
    #=======================Data Export/Import===================================================================
    path('exportData/', exportCSV, name="dataexport"),
    path('import_organism/', import_excel_organism, name="importOrg"),
    path('import_Taxonomy/', import_excel_taxo, name="importTaxo"),
    
    path('import_dictionary/', import_excel_dict, name="importDict"),
]