from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dcell.views import  (Cell_ListView,Cell_CardView,createCell, detailCell, updateCell, Cell_DeleteView,
                     CellBatch_ListView, createCellBatch, CellBatch_DeleteView, CellBatch_UpdateView, 
                     CellBatchStock_ListView, CellcreateStock, CellupdateStock, CellstockList, CellBatchStock_DeleteView,)
from dcell.cell_upload_views import Import_CellView
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
    path('cellbatch_list', CellBatch_ListView.as_view(), name="cellbatch_list"),
    path('createBatch/<str:cell_id>/', createCellBatch, name="cell_batch_create"),
    path('updateBat/<str:pk>', CellBatch_UpdateView.as_view(), name="cell_batch_update"),
    path('deleteBat/<str:pk>', CellBatch_DeleteView.as_view(), name="cell_batch_delete"),

    # CellBatch Stock
    path('stocklist/<str:pk>', CellstockList, name="cell_stock_list"),
    path('stocklist', CellBatchStock_ListView.as_view(), name="cell_stock_list_overview"),
    path('createStock/<str:cellbatch_id>/', CellcreateStock, name="cell_stock_create"),
    path('updateStock/<str:pk>/', CellupdateStock, name="cell_stock_update"),
    path('deleteStock/<str:pk>/', CellBatchStock_DeleteView.as_view(), name="cell_stock_delete"),

    path('import-cell/', Import_CellView.as_view(), name='import-cell'),

    # path('pivottable/<str:pk>', pivottable, name="pivottable"),
  
]