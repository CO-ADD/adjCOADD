from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (exportCSV, import_excel_taxo,import_excel_dict, import_excel_organism, createOrganism, detailOrganism,
                        updateOrganism, updateBatch, TaxonomyCardView,TaxonomyListView,detailTaxonomy, createTaxonomy,updateTaxonomy, deleteTaxonomy, 
                       OrganismListView,OrganismCardView, BatchCardView, createBatch, deleteOrganism, deleteBatch, createStock, updateStock, stockList, deleteStock) #StockListView 
from .utils import search_organism, search_organism_id


urlpatterns = [
    path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
    path('taxonomy_list', TaxonomyListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', detailTaxonomy, name="taxo_detail"),
    path('organism_card', OrganismCardView.as_view(), name="org_card"),
    path('organism_list', OrganismListView.as_view(), name="org_list"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
    path('organism-batch_card', BatchCardView.as_view(), name="batch_card"),
    path('stocklist/<str:pk>', stockList, name="stock_list"),

    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),
    path('createOrg/', createOrganism, name="org_create"),
    path('createTaxo/', createTaxonomy, name="taxo_create"),
    path('createBatch/', createBatch, name="batch_create"),
    path('createStock/', createStock, name="stock_create"),
    

    path('updateTax/<slug:slug>', updateTaxonomy, name="taxonomy_update"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),
    path('updateBat/<str:pk>', updateBatch, name="batch_update"),
    path('organism/updateStock/<str:pk>', updateStock, name="stock_update"),
   

    path('deleteOrg/<str:pk>', deleteOrganism, name="organism_delete"),
    path('deleteTax/<slug:slug>', deleteTaxonomy, name="taxonomy_delete"),
    path('deleteBat/<str:pk>', deleteBatch, name="batch_delete"),
    path('organism/deleteStock/<str:pk>', deleteStock, name="stock_delete"),
  
    #=======================Data Export/Import ===================================================================
    path('exportData/', exportCSV, name="dataexport"),
    path('import_organism/', import_excel_organism, name="importOrg"),
    path('import_Taxonomy/', import_excel_taxo, name="importTaxo"),
    
    path('import_dictionary/', import_excel_dict, name="importDict"),
]