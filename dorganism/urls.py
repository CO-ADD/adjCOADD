from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (TaxonomyCardView,TaxonomyListView,detailTaxonomy,TaxonomyCreateView, TaxonomyUpdateView, deleteTaxonomy, 
                     createOrganism, detailOrganism, updateOrganism, deleteOrganism, OrganismListView,OrganismCardView, 
                     BatchUpdateView, BatchCardView, createBatch, deleteBatch,
                     createStock, updateStock, stockList, deleteStock,
                     updateCulture, createCulture, deleteCulture, pivottable) #StockListView 
from .utils import search_organism, search_organism_id


urlpatterns = [
    # Taxonomy
    path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
    path('taxonomy_list', TaxonomyListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', detailTaxonomy, name="taxo_detail"),
    path('createTaxo/', TaxonomyCreateView.as_view(), name="taxo_create"),
    path('updateTax/<slug:slug>', TaxonomyUpdateView.as_view(), name="taxonomy_update"),
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
    path('createBatch/<str:organism_id>/', createBatch, name="batch_create"),
    path('updateBat/<str:pk>', BatchUpdateView.as_view(), name="batch_update"),
    path('deleteBat/<str:pk>', deleteBatch, name="batch_delete"),

    # OrgBatch Stock
    path('stocklist/<str:pk>', stockList, name="stock_list"),
    path('createStock/<str:orgbatch_id>/', createStock, name="stock_create"),
    path('updateStock/<str:pk>/', updateStock, name="stock_update"),
    path('deleteStock/<str:pk>/', deleteStock, name="stock_delete"),

    # Organism Culture 
    path('createCulture/<str:organism_id>/', createCulture, name="culture_create"),
    path('updateCulture/<str:pk>', updateCulture, name="cultr_update"),
    path('deleteCulture/<str:pk>', deleteCulture, name="cultr_delete"),
  
    #Json search organism
    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),

    path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]