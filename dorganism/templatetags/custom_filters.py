from django import template
import re
from dorganism.models import OrgBatch_Stock
from adjcoadd.constants import LinkList

register = template.Library()

@register.filter
def count_filtered_stock(object_batch):
    return OrgBatch_Stock.objects.filter(orgbatch_id=object_batch, astatus__gte=0, n_left__gt=1).count()

@register.filter
def is_dict(value):
    return isinstance(value, dict)

@register.filter
def is_list(value):
    return isinstance(value, list)

@register.filter
def to_valid_selector(value):
    return "a" + re.sub(r"[\s\.\#\[\]\(\)\+\>\~\=\'\*\^\$]", "_", str(value))

@register.filter
def to_int(value):
    if value != "":
        return int(value)
    else:
        return 0

@register.filter
def get_linkname(value):
    if value in LinkList.keys():
        return LinkList[value]
    else:
        return None

