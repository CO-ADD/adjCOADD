
from django.contrib import admin

from .models import Organisms, Taxonomy
# Register your models here.

class OrganismsAdmin(admin.ModelAdmin):
    list_display=("Organism_ID","Organism_Name")
class TaxonomyAdmin(admin.ModelAdmin):
    pass
admin.site.register(Organisms, OrganismsAdmin)
admin.site.register(Taxonomy, TaxonomyAdmin)

# admin.site.register(Groupfilter)