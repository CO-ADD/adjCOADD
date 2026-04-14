from django.contrib import admin

from .models import  Organisation, Collab_User, Collab_Group, Collab_Membership
# Register your models here.
admin.site.register(Organisation)
admin.site.register(Collab_User)
admin.site.register(Collab_Group)
admin.site.register(Collab_Membership)