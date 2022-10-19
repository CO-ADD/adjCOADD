
from django.contrib import admin
from django.contrib.auth.admin import UserAdmin
from .models import User, ApplicationUser #Dictionaries,
# Register your models here.

class UserAdmin(admin.ModelAdmin):
    list_display=("username","email","is_staff")

# class DictAdmin(admin.ModelAdmin):
    # list_display=("id","Dictionary_ID")
    

admin.site.register(User, UserAdmin)

# admin.site.register(Dictionaries, DictAdmin)
admin.site.register(ApplicationUser)
