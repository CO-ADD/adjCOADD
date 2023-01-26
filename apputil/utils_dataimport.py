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
from apputil.utils import Validation_Log, instance_dict
from ddrug.Vitek import *
from ddrug.models import VITEK_Card, VITEK_ID, VITEK_AST
from .utils import instance_dict, Validation_Log

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
        # if self.max_size is not None and data.size > self.max_size:
        #     params = {
        #         'max_size': filesizeformat(self.max_size), 
        #         'size': filesizeformat(data.size),
        #     }
        #     raise ValidationError(self.error_messages['max_size'],
        #                            'max_size', params)

        # if self.min_size is not None and data.size < self.min_size:
        #     params = {
        #         'min_size': filesizeformat(self.min_size),
        #         'size': filesizeformat(data.size)
        #     }
        #     raise ValidationError(self.error_messages['min_size'], 
        #                            'min_size', params)

        if self.content_types:
            content_type = magic.from_buffer(data.read(), mime=True)
            data.seek(0)

            if content_type not in self.content_types:
                params = { 'content_type': content_type }
                raise ValidationError(self.error_messages['content_type'],
                                   'content_type', params)

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
            object_list.append(Dictionary( dict_value=dbframe.dict_value, dict_class=dbframe.dict_class, dict_desc =dbframe.dict_desc ))
        except Exception as err:
            print(err)
            return err
    print(object_list)
    return object_list

# =========================================================
from django.conf import settings
@transaction.atomic
def import_excel(file_path, data_model):
    print('importing....')
    model=data_model
    excel_file=file_path
    print(excel_file)
    # object_list=[]
    file_name=file_path.split("/")[2]
    dirname=settings.MEDIA_ROOT
    
    if model=='Vitek':
        print("working on Vitek pdf")
        return process_VitekPDF(DirName=dirname, PdfName=file_name) #DirName,PdfName
    
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



# ==========================File Process=====================================

from django.core.files.storage import FileSystemStorage
from django.views import View
from django.views.generic.edit import FormView
from django import forms
import json
from django.core import serializers
import os
# from .utils_dataimport import FileValidator, uploadedfile_process
from django.core.exceptions import ValidationError
from django.db import transaction, IntegrityError
from pathlib import Path
from django.conf import settings
from .utils import instance_dict, Validation_Log, SuperUserRequiredMixin

from ddrug.models import VITEK_Card, VITEK_ID, VITEK_AST


# set filefield Validator
# validate_file = FileValidator(#max_size=1024 * 100, 
#                              content_types=('text/csv', 'application/pdf','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
# create array for files if infected
# infected_files = []
# setup unix socket to scan stream
# cd = clamd.ClamdUnixSocket()

if settings.DEVELOPMENT:
    path='uploads'
else:
    Base_dir = Path(__file__).resolve().parent.parent.parent
    path=os.path.abspath(os.path.join(Base_dir, 'uploads'))

    # #delete task

def delete_file(file_path):
    file_name=file_path.split("/")[2]
    print(file_name)
    file_full_path=os.path.join(settings.MEDIA_ROOT, file_name)
    print(file_full_path)
    try:
        os.unlink(file_full_path)
        print("removed!")
    except Exception as err:
        print(err)


class FileUploadForm(SuperUserRequiredMixin, forms.Form):
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True,}), )#validators=[validate_file])
    

class Importhandler(SuperUserRequiredMixin, FormView):
    
    form_class=FileUploadForm
    file_url=[]
    data_list=[]
    data_model='default'
    success_url="default"
    template_name='default'
 
    # def get(self, request):
    #     form = self.form_class
    #     for f in os.listdir(path):
    #         print(f)
    #     return render(request, 'ddrug/importdata_vitek.html', { 'form': form, })
    
   
    # def form_valid(self, form):
    #     self.data_model=self.request.POST.get('file_data')
    #     myfiles=self.request.FILES.getlist('file_field')
    #     self.file_url=[]
    #     # This method is called when valid form data has been POSTed.
    #     # It should return an HttpResponse.
    #     form.send_email()
    #     return super().form_valid(form)
    
    # def post(self, request):
    #     form = self.form_class(request.POST, request.FILES)
    #     context = {}
    #     context['form'] = form
    #     vLog = Validation_Log('VITEK PDF')
    #     kwargs={}
    #     kwargs['user']=request.user
    #     vCards=[]
    #     vID=[]
    #     vAST=[]
    #     self.data_model=request.POST.get('file_data')
    #     myfiles=request.FILES.getlist('file_field')
    #     self.file_url=[]
    #     try:
    #         if form.is_valid():
    #             print(myfiles)
    #             # scan_results = cd.instream(myfile) # scan_results['stream'][0] == 'OK' or 'FOUND'
    #             for f in myfiles:
    #                 fs=FileSystemStorage()
    #                 filename=fs.save(f.name, f)
    #                 self.file_url.append(fs.url(filename))
    #                 print(self.file_url)
    #                 try:
    #                     vCards,vID,vAst=import_excel(fs.url(filename), self.data_model)
    #                     print("file checked")
           
    #                 except Exception as err:
    #                     print(f"uploaderror is {err}")

    #                 # return messages.warning(request, f'There is {err} error, upload again')
    #             context['file_pathlist']=self.file_url
    #             context['data_model']=self.data_model
    #             # return render(request,'ddrug/importdata_vitek.html', context)
    #         else:
    #             messages.warning(request, f'There is {form.errors} error, upload again')          

    #     except Exception as err:
    #         messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')
        
    #     if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
    #         print(f'self.file_url is : {self.file_url}')
    #         table_name=[","]
    #         validate_result=[","]
    #         file_report=[","]
    #         process_name=request.POST.get('type')
    #         file_pathlist=request.POST.getlist("filepathlist[]")
    #         print(f'selected : {file_pathlist}')
    #         data_model=request.POST.get("datamodel")
    #             # uploadedfile_process(request, table_name,validate_result, file_report, process_name, f, data_model, vLog)
    #         if process_name=='Validation':
    #             for f in file_pathlist:
    #             # models objects list coming from parsed file
    #                 vCards,vID,vAst=import_excel(f, data_model)
    #             # validating each objectslist 
    #                 if vCards:
    #                     table_name.append("Vitek_card")
    #                     for e in vCards:
    #                         djCard=VITEK_Card.check_from_dict(e, vLog)
    #                         validate_result.append(f"CARD status: {djCard.validStatus}")
    #                         file_report.append(str(vLog.show()))
        
    #                 if vID:
    #                     table_name.append("Vitek_id")
    #                     for e in vID:
    #                         djID=VITEK_ID.check_from_dict(e, vLog)
    #                         validate_result.append(f"ID status: {djID.validStatus}")
    #                         file_report.append(str(vLog.show()))
        
    #                 if vAst:
    #                     table_name.append("Vitek_ast")
    #                     for e in vAst:
    #                         djAst=VITEK_AST.check_from_dict(e, vLog)
    #                         validate_result.append(f"AST status: {djAst.validStatus}")
    #                         file_report.append(str(vLog.show()))               
        
    #             return JsonResponse({"table name":"VITEK".join(table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report)})                                   
       
    #         elif process_name=='Cancel':
    #             # Cancel Task
    #             for f in file_pathlist:
    #                 delete_file(file_path=f)
    #             return JsonResponse({"table name":(",").join( table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report)})                                   
    
    #         elif process_name=='DB_Validation':
    #             # import data to DB
    #             for f in file_pathlist:             
    #                 vCards, vID, vAst=import_excel(f, data_model)
    #     # validating each objectslist 
    #                 if vCards:
    #                     table_name.append("Vitek_card")
    #                     for e in vCards:
    #                         djCard=VITEK_Card.check_from_dict(e, vLog)
    #                         if djCard.validStatus:
    #                             try:
    #                                 djCard.save(**kwargs)
    #                             except Exception as err:
    #                                 validate_result.append(f"catch Exception CARD {err}")
    #                         validate_result.append(f"CARD status: {djCard.validStatus}")
    #                         file_report.append(str(vLog.show()))
                
    #                 if vID:
    #                     table_name.append("Vitek_id")
    #                     for e in vID:
    #                         djID=VITEK_ID.check_from_dict(e, vLog)
    #                         if djID.validStatus:
    #                             try:
    #                                 djID.save(**kwargs)
    #                             except Exception as err:
    #                                 validate_result.append(f"catch Exception ID {err}") 
    #                         validate_result.append(f"ID status: {djID.validStatus}")
    #                         file_report.append(str(vLog.show()))
                
    #                 if vAst:
    #                     table_name.append("Vitek_ast")
    #                     for e in vAst:
    #                         djAst=VITEK_AST.check_from_dict(e, vLog)
    #                         if djAst.validStatus:
    #                             try:
    #                                 djAst.save(**kwargs)
    #                             except Exception as err:
    #                                 validate_result.append(f"catch Exception Ast {err}")    
    #                         validate_result.append(f"AST status: {djAst.validStatus}")
    #                         file_report.append(str(vLog.show()))               
                
    #             return JsonResponse({"table name":"VITEK".join(table_name), 'validate_result':(",").join(validate_result), 'file_report':(",").join(file_report), 'status':"Data Saved!"})

           
        # return render(request, 'ddrug/importdata_vitek.html', context)




