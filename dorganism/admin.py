
from django.contrib import admin


from .models import  Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
# Register your models here.
admin.site.register(Organism)
admin.site.register(Organism_Batch)
admin.site.register(OrgBatch_Stock)
admin.site.register(Organism_Culture)
admin.site.register(Taxonomy)


