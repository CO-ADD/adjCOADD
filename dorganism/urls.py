from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dorganism.views import  (Taxonomy_ListView,Taxonomy_CardView,Taxonomy_DetailView,Taxonomy_CreateView, Taxonomy_UpdateView, Taxonomy_DeleteView, 
                     Organism_ListView,Organism_CardView,Organism_CreateView, Organism_DetailView, Organism_UpdateView, Organism_DeleteView,
                     OrgBatch_ListView, OrgBatch_CreateView, OrgBatch_DeleteView, OrgBatch_UpdateView, 
                     OrgBatchStock_ListView, OrgBatchStock_CreateView, OrgBatchStock_UpdateView, OrgBatchStock_DetailView, OrgBatchStock_DeleteView,
                     OrgCulture_UpdateView, OrgCulture_CreateView, OrgCulture_DeleteView,
                     OrgBatchImg_DeleteView,OrgBatchImg_CreateView) 
from dorganism.utils.utils import search_organism, search_organism_id


urlpatterns = [
    # Taxonomy
    path('taxonomy_card', Taxonomy_CardView.as_view(), name="taxo_card"),
    path('taxonomy_list', Taxonomy_ListView.as_view(), name="taxo_list"),
    path('taxonomy/<slug:slug>', Taxonomy_DetailView, name="taxo_detail"),
    path('createTaxo/', Taxonomy_CreateView.as_view(), name="taxo_create"),
    path('updateTax/<slug:slug>', Taxonomy_UpdateView.as_view(), name="taxonomy_update"),
    path('deleteTax/<slug:slug>', Taxonomy_DeleteView.as_view(), name="taxonomy_delete"),

    # Organism 
    path('organism_card', Organism_CardView.as_view(), name="org_card"),
    path('organism_list', Organism_ListView.as_view(), name="org_list"),
    path('organism/<str:pk>', Organism_DetailView, name="org_detail"),
    path('createOrg/', Organism_CreateView, name="org_create"),
    path('updateOrg/<str:pk>', Organism_UpdateView, name="organism_update"),
    path('deleteOrg/<str:pk>', Organism_DeleteView.as_view(), name="organism_delete"),

    # OrgBatch
    path('batchlist', OrgBatch_ListView.as_view(), name="org_batch_list"),
    path('createBatch/<str:organism_id>/', OrgBatch_CreateView, name="org_batch_create"),
    path('updateBat/<str:pk>', OrgBatch_UpdateView.as_view(), name="org_batch_update"),
    path('deleteBat/<str:pk>', OrgBatch_DeleteView.as_view(), name="org_batch_delete"),
    
    # OrgBatch Images
    path('createbatchimg/<str:pk>', OrgBatchImg_CreateView.as_view(), name="batchimg_create"),
    path('deleteBatimg/<str:pk>', OrgBatchImg_DeleteView.as_view(), name="batchimg_delete"),

    # OrgBatch Stock
    path('stocklist/<str:pk>', OrgBatchStock_DetailView, name="org_stock_list"),
    path('stocklist', OrgBatchStock_ListView.as_view(), name="org_stock_list_overview"),
    path('createStock/<str:orgbatch_id>/', OrgBatchStock_CreateView, name="org_stock_create"),
    path('updateStock/<str:pk>/', OrgBatchStock_UpdateView, name="org_stock_update"),
    path('deleteStock/<str:pk>/', OrgBatchStock_DeleteView.as_view(), name="org_stock_delete"),

    # Organism Culture 
    path('createCulture/<str:organism_id>/', OrgCulture_CreateView, name="culture_create"),
    path('updateCulture/<str:pk>', OrgCulture_UpdateView.as_view(), name="culture_update"),
    path('deleteCulture/<str:pk>', OrgCulture_DeleteView.as_view(), name="culture_delete"),
  
    #Json search organism
    path('search_organism/', search_organism, name="search_organism"),
    path('search_organism_id/', search_organism_id, name="search_organism_id"),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]