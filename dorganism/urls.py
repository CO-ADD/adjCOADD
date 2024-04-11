from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dorganism.views import  (Taxonomy_ListView,Taxonomy_CardView,detailTaxonomy,Taxonomy_CreateView, Taxonomy_UpdateView, Taxonomy_DeleteView, 
                     Organism_ListView,Organism_CardView,createOrganism, detailOrganism, updateOrganism, Organism_DeleteView,
                     OrgBatch_ListView, createBatch, OrgBatch_DeleteView, OrgBatch_UpdateView, 
                     OrgBatchStock_ListView, createStock, updateStock, stockList, OrgBatchStock_DeleteView,
                     OrgCulture_UpdateView, createCulture, OrgCulture_DeleteView,
                     OrgBatchImg_DeleteView,OrgBatchImg_CreateView) 
from dorganism.utils.utils import search_organism, search_organism_id


urlpatterns = [
    # Taxonomy
    path('taxonomy_card', Taxonomy_CardView.as_view(), name="taxo_card"),
    path('taxonomy_list', Taxonomy_ListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', detailTaxonomy, name="taxo_detail"),
    path('createTaxo/', Taxonomy_CreateView.as_view(), name="taxo_create"),
    path('updateTax/<slug:slug>', Taxonomy_UpdateView.as_view(), name="taxonomy_update"),
    path('deleteTax/<slug:slug>', Taxonomy_DeleteView.as_view(), name="taxonomy_delete"),

    # Organism 
    path('organism_card', Organism_CardView.as_view(), name="org_card"),
    path('organism_list', Organism_ListView.as_view(), name="org_list"),
    path('organism/<str:pk>', detailOrganism, name="org_detail"),
    path('createOrg/', createOrganism, name="org_create"),
    path('updateOrg/<str:pk>', updateOrganism, name="organism_update"),
    path('deleteOrg/<str:pk>', Organism_DeleteView.as_view(), name="organism_delete"),

    # OrgBatch
    path('batchlist', OrgBatch_ListView.as_view(), name="batch_list"),
    path('createBatch/<str:organism_id>/', createBatch, name="org_batch_create"),
    path('updateBat/<str:pk>', OrgBatch_UpdateView.as_view(), name="batch_update"),
    path('deleteBat/<str:pk>', OrgBatch_DeleteView.as_view(), name="batch_delete"),
    path('createbatchimg/<str:pk>', OrgBatchImg_CreateView.as_view(), name="batchimg_create"),
    path('deleteBatimg/<str:pk>', OrgBatchImg_DeleteView.as_view(), name="batchimg_delete"),

    # OrgBatch Stock
    path('stocklist/<str:pk>', stockList, name="stock_list"),
    path('stocklist', OrgBatchStock_ListView.as_view(), name="org_stock_list_overview"),
    path('createStock/<str:orgbatch_id>/', createStock, name="stock_create"),
    path('updateStock/<str:pk>/', updateStock, name="stock_update"),
    path('deleteStock/<str:pk>/', OrgBatchStock_DeleteView.as_view(), name="stock_delete"),

    # Organism Culture 
    path('createCulture/<str:organism_id>/', createCulture, name="culture_create"),
    path('updateCulture/<str:pk>', OrgCulture_UpdateView.as_view(), name="culture_update"),
    path('deleteCulture/<str:pk>', OrgCulture_DeleteView.as_view(), name="culture_delete"),
  
    #Json search organism
    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]