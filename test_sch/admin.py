from django.contrib import admin

from .models import Taxo
# Register your models here.


class TaxoAdmin(admin.ModelAdmin):
    pass
admin.site.register(Taxo, TaxoAdmin)
