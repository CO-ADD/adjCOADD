from django import template
import re
from dorganism.models import OrgBatch_Stock
from dcell.models import CellBatch_Stock
from adjcoadd.constants import LinkList

register = template.Library()

@register.filter
def count_filtered_org_stock(object_batch):
    return OrgBatch_Stock.objects.filter(orgbatch_id=object_batch, astatus__gte=0, n_left__gt=0).count()

@register.filter
def count_filtered_cell_stock(object_batch):
    return CellBatch_Stock.objects.filter(cellbatch_id=object_batch, astatus__gte=0, n_left__gt=0).count()

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
def get_linkname(value, arg1=None, arg2=None):
    if value in LinkList.keys():        
        arg1_str = str(arg1) if arg1 else ''
        arg2_str = str(arg1) if arg2 else ''
        return LinkList[value].replace('{VALUE1}', arg1_str).replace('{VALUE2}', arg2_str)
    else:
        return None

