
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import  ApplicationUser, Dictionary, ApplicationLog
# Register your models here.

class UserAdmin(admin.ModelAdmin):
    list_display=("username","email","is_staff", "name")

class DictAdmin(admin.ModelAdmin):
    def save(self, request, *args, **kwargs):
        kwargs['user']=request.user
        super().save(*args, **kwargs)
    

admin.site.register(ApplicationUser, UserAdmin)

admin.site.register(Dictionary, DictAdmin)

admin.site.register(ApplicationLog)
