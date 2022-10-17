
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import User, Groupfilter, ApplicationDictionary,ApplicationUser
# Register your models here.

class UserAdmin(admin.ModelAdmin):
    list_display=("username","email","is_staff")

admin.site.register(User, UserAdmin)
admin.site.register(Groupfilter)
admin.site.register(ApplicationDictionary)
admin.site.register(ApplicationUser)
