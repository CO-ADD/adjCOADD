import os
import re
import unicodedata
import django_filters
import pandas as pd

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction

from dorganism.models import Taxonomy
from apputil.models import Dictionary

# -----------------------Start Utility Functions-----------------------------------

def has_special_char(text: str) -> bool:
    return any(c for c in text if not c.isalnum() and not c.isspace())
# ----------------------------------------------------------------------------------

@transaction.atomic
def import_excel(file_path, data_model):
    print('importing....')
    model=data_model
    if model == 'Taxonomy':
        # Validation Taxonomy data import with model 
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
    else:
        try:
            raise Exception('No Model Found, Please Choose Model: Taxonomy, Organism...')
        except Exception  as err:
            return err  
    # return 'something wrong'

