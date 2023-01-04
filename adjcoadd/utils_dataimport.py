import os
import re
import unicodedata
import django_filters
import pandas as pd

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture
from apputil.models import Dictionary

# -----------------------Start Utility Functions-----------------------------------


# ==============Uploading File Validators==============================
# import magic

from django.utils.deconstruct import deconstructible
from django.template.defaultfilters import filesizeformat
from django.core.exceptions import ValidationError


@deconstructible
class FileValidator(object):
    error_messages = {
     'max_size': ("Ensure this file size is not greater than %(max_size)s."
                  " Your file size is %(size)s."),
     'min_size': ("Ensure this file size is not less than %(min_size)s. "
                  "Your file size is %(size)s."),
     'content_type': "Files of type %(content_type)s are not supported.",
    }

    def __init__(self, max_size=None, min_size=None, content_types=()):
        self.max_size = max_size
        self.min_size = min_size
        self.content_types = content_types

    def __call__(self, data):
        if self.max_size is not None and data.size > self.max_size:
            params = {
                'max_size': filesizeformat(self.max_size), 
                'size': filesizeformat(data.size),
            }
            raise ValidationError(self.error_messages['max_size'],
                                   'max_size', params)

        if self.min_size is not None and data.size < self.min_size:
            params = {
                'min_size': filesizeformat(self.min_size),
                'size': filesizeformat(data.size)
            }
            raise ValidationError(self.error_messages['min_size'], 
                                   'min_size', params)

        # if self.content_types:
        #     content_type = magic.from_buffer(data.read(), mime=True)
        #     data.seek(0)

        #     if content_type not in self.content_types:
        #         params = { 'content_type': content_type }
        #         raise ValidationError(self.error_messages['content_type'],
        #                            'content_type', params)

    def __eq__(self, other):
        return (
            isinstance(other, FileValidator) and
            self.max_size == other.max_size and
            self.min_size == other.min_size and
            self.content_types == other.content_types
        )


def has_special_char(text: str) -> bool:
    return any(c for c in text if not c.isalnum() and not c.isspace())
# ----------------------------------------------------------------------------------
# =========================================================
def validate_Taxonomy(dbframe):
    object_list=[]
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
                    object_list.append(Taxonomy(organism_name=dbframe.ORGANISM_NAME, other_names=dbframe.ORGANISM_NAME_OTHER, code=dbframe.ORGANISM_CODE, 
                        org_class=class_fkey, tax_id=dbframe.TAX_ID, parent_tax_id=dbframe.PARENT_TAX_ID, 
                        tax_rank=dbframe.TAX_RANK, division=division_fkey, lineage=linea, 
                        ))
                except Exception as err:
                    print(err)
                    # obj.save()
            
            except Exception as err:
                print(err)
                return err
    return object_list

# =========================================================
def validate_Organism(dbframe):
    object_list=[]
    for dbframe in dbframe.itertuples():
        try:
            taxID=int('0'+dbframe[22])
            screen_panel=dbframe[26].split(';')
            organism_fkey=Taxonomy.objects.filter(organism_name=dbframe[1])
            print(organism_fkey[0])   
            try:
                object_list.append(Organism(organism_id=dbframe[0], organism_name=organism_fkey[0],  strain_id=dbframe[3], 
                                    strain_code=dbframe[5], strain_notes=dbframe[7], 
                                    strain_tissue=dbframe[25], strain_type=dbframe[4], sequence=dbframe[28], sequence_link=dbframe[29], 
                                    strain_panel=screen_panel, 
                                    tax_id =taxID,risk_group=dbframe[9], pathogen_group =dbframe[10],import_permit =dbframe[12],bio_approval =dbframe[23],special_precaution =dbframe[24],lab_restriction =dbframe[27],mta_document =dbframe[31],
                                    mta_status =dbframe[32],oxygen_pref =dbframe[13],atmosphere_pref ='containSpecialCHA', nutrient_pref =dbframe[15],biofilm_pref =dbframe[16],))
            except Exception as err:
                print(err)
        except Exception as err:
            print(err)
            return err
    return object_list
# =========================================================
def validate_Dictionary(dbframe):
    object_list=[]
    for dbframe in dbframe.itertuples():
        try:
            object_list.append(Dictionary(dict_class=dbframe.Class, dict_value=dbframe.Term, dict_desc =dbframe.Name, ))
        except Exception as err:
            print(err)
            return err
    return object_list

# =========================================================
@transaction.atomic
def import_excel(file_path, data_model):
    print('importing....')
    model=data_model
    excel_file=file_path
    print(excel_file)
    # object_list=[]
    
    try:
        exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
    except Exception as err:
        return err
    dbframe=exmpexceldata

    if model == 'Taxonomy':
        return validate_Taxonomy(dbframe)
    elif model == 'Organism':
        return validate_Organism(dbframe)
    elif model=='Dictionary':
        return validate_Dictionary(dbframe)

    else:
        raise Exception('No Model Found, Please Choose Model: Taxonomy, Organism...')
        return err  
    # return 'something wrong'