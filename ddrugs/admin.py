from django.contrib import admin

from .models import Drugbank
# Register your models here.


class DrugbankAdmin(admin.ModelAdmin):
    pass
admin.site.register(Drugbank, DrugbankAdmin)
