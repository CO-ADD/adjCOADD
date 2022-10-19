
from django.contrib import admin

from .models import Organisms, Taxonomy
# Register your models here.

class OrganismsAdmin(admin.ModelAdmin):
    pass
class TaxonomyAdmin(admin.ModelAdmin):
    pass




admin.site.register(Organisms, OrganismsAdmin)
admin.site.register(Taxonomy, TaxonomyAdmin)


# # admin.site.register(Groupfilter)
# from django.contrib.admin.models import LogEntry

# # @admin.register(LogEntry)
# class LogEntryAdmin(admin.ModelAdmin):
#     # to have a date-based drilldown navigation in the admin page
#     date_hierarchy = 'action_time'

#     # to filter the resultes by users, content types and action flags
#     list_filter = [
#         'user',
#         'content_type',
#         'action_flag'
#     ]

#     # when searching the user will be able to search in both object_repr and change_message
#     search_fields = [
#         'object_repr',
#         'change_message'
#     ]

#     list_display = [
#         'action_time',
#         'user',
#         'content_type',
#         'action_flag',
#     ]