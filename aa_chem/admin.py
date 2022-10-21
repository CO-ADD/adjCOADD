
from django.contrib import admin

from .models import  Mytest,Organisms, Taxonomy
# Register your models here.

# class OrganismsAdmin(admin.ModelAdmin):    

#     def formfield_for_choice_field(self, db_field, request, **kwargs):
#         if db_field.name == "Strain_Type":
#             kwargs['choices'] = (
#                 ('accepted', 'Accepted'),
#                 ('denied', 'Denied'),
#             )
#             if request.user.is_superuser:
#                 kwargs['choices'] += (('ready', 'Ready for deployment'),)
#         return super(OrganismsAdmin, self).formfield_for_choice_field(db_field, request, **kwargs)



# class TaxonomyAdmin(admin.ModelAdmin):
#     pass



    

admin.site.register(Organisms)
admin.site.register(Taxonomy)
admin.site.register(Mytest)

