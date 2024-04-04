from django.contrib import admin

from .models import  Cell, Cell_Batch, CellBatch_Stock
# Register your models here.
admin.site.register(Cell)
admin.site.register(Cell_Batch)
admin.site.register(CellBatch_Stock)

