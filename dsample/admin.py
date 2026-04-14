from django.contrib import admin

# Register your models here.
from .models import  Project, Project_Membership,Library,Compound_Batch,COADD_Compound,ABase_Compound,Library_Compound
# Register your models here.
admin.site.register(Project)
admin.site.register(Project_Membership)
admin.site.register(Library)
admin.site.register(Compound_Batch)
admin.site.register(COADD_Compound)
admin.site.register(ABase_Compound)
admin.site.register(Library_Compound)