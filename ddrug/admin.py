from django.contrib import admin
from .models import Drug, VITEK_Card, VITEK_AST,VITEK_ID
# Register your models here.

class DrugAdmin(admin.ModelAdmin):
    list_display=("drug_id", "drug_name", "drug_othernames", "acreated")


admin.site.register(Drug, DrugAdmin)

admin.site.register(VITEK_Card)

admin.site.register(VITEK_AST)

admin.site.register(VITEK_ID)


