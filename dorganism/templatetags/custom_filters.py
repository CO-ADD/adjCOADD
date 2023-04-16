from django import template
from dorganism.models import OrgBatch_Stock

register = template.Library()

@register.filter
def count_filtered_instances(object_batch):
    return OrgBatch_Stock.objects.filter(orgbatch_id=object_batch, astatus__gte=0, n_left__gte=1).count()