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

# =========================================================
# @login_required
# @user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
# def import_excel_dict(req):
#     print('importing....')
#     try:
#         if req.method=='POST' and req.FILES['myfile']:
#             myfile=req.FILES['myfile']
#             fs=FileSystemStorage()
#             filename=fs.save(myfile.name, myfile)
#             uploaded_file_url=fs.url(filename)
#             excel_file=uploaded_file_url
#             print(excel_file)
#             exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
#             print(type(exmpexceldata))
#             dbframe=exmpexceldata
#             for dbframe in dbframe.itertuples():                   
#                 obj, created=Dictionary.objects.get_or_create(dict_class=dbframe.Class, dict_value=dbframe.Term, dict_desc =dbframe.Name, )
#                 print(type(obj))
          
#             return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {'uploaded_file_url': uploaded_file_url})
#     except Exception as err:
#         print(f'import failed because {err}')
#     return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {})

# #==================================================================import Organism================================================
# @login_required
# @user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
# def import_excel_organism(req):
#     print('importing....')
#     try:
#         if req.method=='POST' and req.FILES['myfile']:
#             myfile=req.FILES['myfile']
#             fs=FileSystemStorage()
#             filename=fs.save(myfile.name, myfile)
#             uploaded_file_url=fs.url(filename)
#             excel_file=uploaded_file_url
#             print(excel_file)
#             exmpexceldata=pd.read_csv("."+excel_file, )
#             print(exmpexceldata.itertuples)
#             dbframe=exmpexceldata
#             for dbframe in dbframe.itertuples():
#                 taxID=int('0'+dbframe[22])
#                 screen_panel=dbframe[26].split(';')
#                 organism_fkey=Taxonomy.objects.filter(organism_name=dbframe[1])
#                 print(organism_fkey[0])   
#                 try:
#                     obj, created=Organism.objects.get_or_create(organism_id=dbframe[0], organism_name=organism_fkey[0],  strain_id=dbframe[3], 
#                                     strain_code=dbframe[5], strain_notes=dbframe[7], 
#                                     strain_tissue=dbframe[25], strain_type=dbframe[4], sequence=dbframe[28], sequence_link=dbframe[29], 
#                                     strain_panel=screen_panel, 
#                                     tax_id =taxID,risk_group=dbframe[9], pathogen_group =dbframe[10],import_permit =dbframe[12],bio_approval =dbframe[23],special_precaution =dbframe[24],lab_restriction =dbframe[27],mta_document =dbframe[31],
#                                     mta_status =dbframe[32],oxygen_pref =dbframe[13],atmosphere_pref ='containSpecialCHA', nutrient_pref =dbframe[15],biofilm_pref =dbframe[16], acreated_by=req.user )
#                 except Exception as err:
#                     print(err)
#                 # obj.save()
            
#             return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {'uploaded_file_url': uploaded_file_url})
#     except Exception as err:
#         print(err)
#     return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {})