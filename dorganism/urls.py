from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (TaxonomyCardView,TaxonomyListView,detailTaxonomy, createTaxonomy, updateTaxonomy, deleteTaxonomy, 
                     createOrganism, detailOrganism, updateOrganism, deleteOrganism, OrganismListView,OrganismCardView, 
                     updateBatch, BatchCardView, createBatch, deleteBatch,
                     createStock, updateStock, stockList, deleteStock,
                     updateCulture, createCulture, deleteCulture) #StockListView 
from .utils import search_organism, search_organism_id


urlpatterns = [
    # Taxonomy
    path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
    path('taxonomy_list', TaxonomyListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', detailTaxonomy, name="taxo_detail"),
    path('createTaxo/', createTaxonomy, name="taxo_create"),
    path('updateTax/<slug:slug>', updateTaxonomy, name="taxonomy_update"),
    path('deleteTax/<slug:slug>', deleteTaxonomy, name="taxonomy_delete"),

    # Organism 
    path('organism_card', OrganismCardView.as_view(), name="org_card"),
   
    # ------------------------------------------------------------------
    path('organism_list', OrganismListView.as_view(), name="org_list"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
    path('createOrg/', createOrganism, name="org_create"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),
    path('deleteOrg/<str:pk>', deleteOrganism, name="organism_delete"),

    # OrgBatch
    path('organism-batch_card', BatchCardView.as_view(), name="batch_card"),
    path('createBatch/', createBatch, name="batch_create"),
    path('updateBat/<str:pk>', updateBatch, name="batch_update"),
    path('deleteBat/<str:pk>', deleteBatch, name="batch_delete"),

    # OrgBatch Stock
    path('stocklist/<str:pk>', stockList, name="stock_list"),
    path('createStock/', createStock, name="stock_create"),
    path('organism/updateStock/<str:pk>', updateStock, name="stock_update"),
    path('organism/deleteStock/<str:pk>', deleteStock, name="stock_delete"),

    # Organism Culture 
    path('createCulture/', createCulture, name="cultr_create"),
    path('updateCulture/<str:pk>', updateCulture, name="cultr_update"),
    path('deleteCulture/<str:pk>', deleteCulture, name="cultr_delete"),
  
    #Json search organism
    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),
  
]