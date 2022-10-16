
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import User, Groupfilter, ApplicationDictionary
# Register your models here.

class UserAdmin(admin.ModelAdmin):
    list_display=("username","email","is_staff","department")

admin.site.register(User, UserAdmin)
admin.site.register(Groupfilter)
admin.site.register(ApplicationDictionary)
