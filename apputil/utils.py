import os
import re
import unicodedata
import django_filters
import pandas as pd

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction

from .models import Dictionary


# ===================================Dictionary query convert to choice Tuples========================================================================#

def get_DictonaryChoices_byDictClass(ModelName, DictClass, sep='|'):
    options=ModelName.objects.filter(dict_class=DictClass).values('dict_value', 'dict_desc')
    if options:
        choices_values=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]+sep+a[1]) for a in choices_values)
    else:
        choices=(('--', 'empty'),)
    return choices

#-----------------------------------------------------------------------------------
def slugify(value, lower=False, allow_unicode=False):
    """
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize("NFKC", value)
    else:
        value = (
            unicodedata.normalize("NFKD", value)
            .encode("ascii", "ignore")
            .decode("ascii")
        )
    value = re.sub(r"[^\w\s-]", "", value)
    value = re.sub(r"[-\s]+", "-", value).strip("-_")

    if lower:
        return value.lower()
    else:
        return value
#-----------------------------------------------------------------------------------
