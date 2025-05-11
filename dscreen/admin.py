from django.contrib import admin

from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50
# Register your models here.
admin.site.register(Screen_Run)
admin.site.register(Assay)
admin.site.register(AssayData_MIC)
admin.site.register(AssayData_CC50)
admin.site.register(AssayData_HC50)