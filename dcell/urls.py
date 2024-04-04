from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dcell.views import  (Cell_ListView,Cell_CardView,createCell, detailCell, updateCell, Cell_DeleteView,
                     CellBatch_ListView, createBatch, CellBatch_DeleteView, CellBatch_UpdateView, 
                     CellBatchStock_ListView, createStock, updateStock, stockList, CellBatchStock_DeleteView,)
from dorganism.utils.utils import search_organism, search_organism_id

urlpatterns = [

    # Cell
    path('cell_card', Cell_CardView.as_view(), name="cell_card"),
    path('cell_list', Cell_ListView.as_view(), name="cell_list"),
    path('cell/<str:pk>', detailCell, name="cell_detail"),
    path('createCell/', createCell, name="cell_create"),
    path('updateCell/<str:pk>', updateCell, name="cell_update"),
    path('deleteCell/<str:pk>', Cell_DeleteView.as_view(), name="cell_delete"),

    # CellBatch
    path('batchlist', CellBatch_ListView.as_view(), name="batch_list"),
    path('createBatch/<str:organism_id>/', createBatch, name="batch_create"),
    path('updateBat/<str:pk>', CellBatch_UpdateView.as_view(), name="batch_update"),
    path('deleteBat/<str:pk>', CellBatch_DeleteView.as_view(), name="batch_delete"),

    # CellBatch Stock
    path('stocklist/<str:pk>', stockList, name="stock_list"),
    path('stocklist', CellBatchStock_ListView.as_view(), name="stock_list_overview"),
    path('createStock/<str:orgbatch_id>/', createStock, name="stock_create"),
    path('updateStock/<str:pk>/', updateStock, name="stock_update"),
    path('deleteStock/<str:pk>/', CellBatchStock_DeleteView.as_view(), name="stock_delete"),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]