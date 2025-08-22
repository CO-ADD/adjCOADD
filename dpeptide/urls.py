from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path


from dpeptide.views import  (Peptide_ListView,Peptide_CreateView, Peptide_DetailView, Peptide_UpdateView, Peptide_RemoveView,
                    #  PepBatch_ListView, PepBatch_CreateView, PepBatch_RemoveView, PepBatch_UpdateView, 
                    #  CellBatchStock_ListView, CellBatchStock_CreateView, CellBatchStock_UpdateView, CellBatchStock_DetailView, CellBatchStock_RemoveView,
                    #  Cell_Upload_HandlerView
                     )


urlpatterns = [

    # Peptide
    #path('peptide_card', Peptide_CardView.as_view(), name="peptide_card"),
    path('peptide_list', Peptide_ListView.as_view(), name="peptide_list"),
    path('peptide/<str:pk>', Peptide_DetailView, name="peptide_detail"),
    path('createPeptide/', Peptide_CreateView, name="peptide_create"),
    path('updatePeptide/<str:pk>', Peptide_UpdateView, name="peptide_update"),
    path('deletePeptide/<str:pk>', Peptide_RemoveView.as_view(), name="peptide_delete"),

    # PepBatch
    # path('pepbatch_list', PepBatch_ListView.as_view(), name="pepbatch_list"),
    # path('createBatch/<str:pep_id>/', PepBatch_CreateView, name="pep_batch_create"),
    # path('updateBat/<str:pk>', PepBatch_UpdateView.as_view(), name="pep_batch_update"),
    # path('deleteBat/<str:pk>', PepBatch_RemoveView.as_view(), name="pep_batch_delete"),

    # PepBatch Stock
    # path('stocklist/<str:pk>', PepBatchStock_DetailView, name="pep_stock_list"),
    # path('stocklist', PepBatchStock_ListView.as_view(), name="pep_stock_list_overview"),
    # path('createStock/<str:pepbatch_id>/', PepBatchStock_CreateView, name="pep_stock_create"),
    # path('updateStock/<str:pk>/', PepBatchStock_UpdateView, name="pep_stock_update"),
    # path('deleteStock/<str:pk>/', PepBatchStock_RemoveView.as_view(), name="pep_stock_delete"),

    # path('import-pep/', Pep_Upload_HandlerView.as_view(), name='import-pep'),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]