from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dcell.views import  (Cell_ListView,Cell_CardView,Cell_CreateView, Cell_DetailView, Cell_UpdateView, Cell_DeleteView,
                     CellBatch_ListView, CellBatch_CreateView, CellBatch_DeleteView, CellBatch_UpdateView, 
                     CellBatchStock_ListView, CellBatchStock_CreateView, CellBatchStock_UpdateView, CellBatchStock_DetailView, CellBatchStock_DeleteView,
                     Cell_Upload_HandlerView)
from dorganism.utils.utils import search_organism, search_organism_id

urlpatterns = [

    # Cell
    path('cell_card', Cell_CardView.as_view(), name="cell_card"),
    path('cell_list', Cell_ListView.as_view(), name="cell_list"),
    path('cell/<str:pk>', Cell_DetailView, name="cell_detail"),
    path('createCell/', Cell_CreateView, name="cell_create"),
    path('updateCell/<str:pk>', Cell_UpdateView, name="cell_update"),
    path('deleteCell/<str:pk>', Cell_DeleteView.as_view(), name="cell_delete"),

    # CellBatch
    path('cellbatch_list', CellBatch_ListView.as_view(), name="cellbatch_list"),
    path('createBatch/<str:cell_id>/', CellBatch_CreateView, name="cell_batch_create"),
    path('updateBat/<str:pk>', CellBatch_UpdateView.as_view(), name="cell_batch_update"),
    path('deleteBat/<str:pk>', CellBatch_DeleteView.as_view(), name="cell_batch_delete"),

    # CellBatch Stock
    path('stocklist/<str:pk>', CellBatchStock_DetailView, name="cell_stock_list"),
    path('stocklist', CellBatchStock_ListView.as_view(), name="cell_stock_list_overview"),
    path('createStock/<str:cellbatch_id>/', CellBatchStock_CreateView, name="cell_stock_create"),
    path('updateStock/<str:pk>/', CellBatchStock_UpdateView, name="cell_stock_update"),
    path('deleteStock/<str:pk>/', CellBatchStock_DeleteView.as_view(), name="cell_stock_delete"),

    path('import-cell/', Cell_Upload_HandlerView.as_view(), name='import-cell'),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]