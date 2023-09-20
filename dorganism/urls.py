from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (TaxonomyCardView,TaxonomyListView,detailTaxonomy,TaxonomyCreateView, TaxonomyUpdateView, TaxonomyDeleteView, 
                     createOrganism, detailOrganism, updateOrganism, OrganismDeleteView, OrganismListView,OrganismCardView, 
                     BatchUpdateView, createBatch, BatchDeleteView, BatchListView, OrgbatchimgDeleteView,
                     createStock, updateStock, stockList, StockDeleteView,
                     CultureUpdateView, createCulture, CultureDeleteView,
                     OrgBatchStock_ListView, OrgbatchimgCreateView) 
from .utils.utils import search_organism, search_organism_id


urlpatterns = [
    # Taxonomy
    path('taxonomy_card', TaxonomyCardView.as_view(), name="taxo_card"),
    path('taxonomy_list', TaxonomyListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', detailTaxonomy, name="taxo_detail"),
    path('createTaxo/', TaxonomyCreateView.as_view(), name="taxo_create"),
    path('updateTax/<slug:slug>', TaxonomyUpdateView.as_view(), name="taxonomy_update"),
    path('deleteTax/<slug:slug>', TaxonomyDeleteView.as_view(), name="taxonomy_delete"),

    # Organism 
    path('organism_card', OrganismCardView.as_view(), name="org_card"),
   
    # ------------------------------------------------------------------
    path('organism_list', OrganismListView.as_view(), name="org_list"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
    path('createOrg/', createOrganism, name="org_create"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),
    path('deleteOrg/<str:pk>', OrganismDeleteView.as_view(), name="organism_delete"),

    # OrgBatch
    path('batchlist', BatchListView.as_view(), name="batch_list"),
    path('createBatch/<str:organism_id>/', createBatch, name="batch_create"),
    path('updateBat/<str:pk>', BatchUpdateView.as_view(), name="batch_update"),
    path('deleteBat/<str:pk>', BatchDeleteView.as_view(), name="batch_delete"),
    path('createbatchimg/<str:pk>', OrgbatchimgCreateView.as_view(), name="batchimg_create"),
    path('deleteBatimg/<str:pk>', OrgbatchimgDeleteView.as_view(), name="batchimg_delete"),

    # OrgBatch Stock
    path('stocklist/<str:pk>', stockList, name="stock_list"),
    path('stocklist', OrgBatchStock_ListView.as_view(), name="stock_list_overview"),
    path('createStock/<str:orgbatch_id>/', createStock, name="stock_create"),
    path('updateStock/<str:pk>/', updateStock, name="stock_update"),
    path('deleteStock/<str:pk>/', StockDeleteView.as_view(), name="stock_delete"),

    # Organism Culture 
    path('createCulture/<str:organism_id>/', createCulture, name="culture_create"),
    path('updateCulture/<str:pk>', CultureUpdateView.as_view(), name="culture_update"),
    path('deleteCulture/<str:pk>', CultureDeleteView.as_view(), name="culture_delete"),
  
    #Json search organism
    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]