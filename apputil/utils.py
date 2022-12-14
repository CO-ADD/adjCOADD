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
from dorganism.models import *


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
def has_special_char(text: str) -> bool:
    return any(c for c in text if not c.isalnum() and not c.isspace())
# ----------------------------------------------------------------------------------
@transaction.atomic
def import_excel(file_path):
    print('importing....')
    Taxonomy_object_list=[]
    excel_file=file_path
    print(excel_file)
    try:
        exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
    except Exception as err:
        return err
    dbframe=exmpexceldata
    for dbframe in dbframe.itertuples():
        try:
            print(str(dbframe.ORGANISM_CLASS))
            if has_special_char(str(dbframe.ORGANISM_NAME)):
                print("speical")
                raise Exception('datatype')
            class_fkey=Dictionary.objects.filter(dict_value=dbframe.ORGANISM_CLASS)
            if class_fkey:
                class_fkey=class_fkey[0]
            else:
                class_fkey=None
            print(class_fkey)
            division_fkey=Dictionary.objects.filter(dict_value=dbframe.DIVISION)
            if division_fkey:
                division_fkey=division_fkey[0]
            else:
                division_fkey=None
            linea=str(dbframe.LINEAGE).split(";")
            print(division_fkey)
            try:
                Taxonomy_object_list.append(Taxonomy(organism_name=dbframe.ORGANISM_NAME, other_names=dbframe.ORGANISM_NAME_OTHER, code=dbframe.ORGANISM_CODE, 
                    org_class=class_fkey, tax_id=dbframe.TAX_ID, parent_tax_id=dbframe.PARENT_TAX_ID, 
                    tax_rank=dbframe.TAX_RANK, division=division_fkey, lineage=linea, 
                    ))
            except Exception as err:
                print(err)
                # obj.save()
            
        except Exception as err:
            print(err)
            return err
    return Taxonomy_object_list
    # return 'something wrong'

