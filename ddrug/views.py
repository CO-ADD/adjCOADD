import json
import os
from rdkit import Chem
from django_rdkit.models import *
from django_rdkit.config import config
from django_filters.views import FilterView
import pandas as pd
import numpy as np
from django.core.serializers.json import DjangoJSONEncoder
from time import localtime, strftime
import psycopg2
from rest_framework import generics

import logging
logger = logging.getLogger("django")
from django.utils.decorators import classonlymethod
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction, IntegrityError
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.conf import settings
from django.views.generic import ListView, TemplateView
from django.utils.functional import SimpleLazyObject


from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import FilteredListView
from apputil.utils.files_upload import Importhandler, file_location, OverwriteStorage
from apputil.utils.api_filterclass import API_FilteredListView
from apputil.utils.validation_log import Validation_Log
from apputil.utils.views_base import permission_not_granted, SimplecreateView, SimpleupdateView
from adjcoadd.constants import *
from dorganism.models import Organism_Batch
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub
from .utils.molecules import molecule_to_svg, clearIMGfolder, get_mfp2_neighbors
from .utils.vitek import *
from .forms import Drug_form, Drug_filter, Vitekcard_filter, Vitekast_filter, MIC_COADDfilter, MIC_Pubfilter
from .serializers import Drug_Serializer, VITEK_ASTSerializer

# ===================================================================
@login_required   
def smartsQuery(req, pk):
    '''
    MAKE SUBSTRUCTURE QUERY 
    '''
    context={}
    object_=get_object_or_404(Drug, drug_id=pk)
    context["object"]=object_
    # get mol block for an object
    try:
        context["object_mol"]=Chem.MolToMolBlock(object_.smol)
    # convert object to JMSE regonized form
        m="\\n".join(context["object_mol"].split("\n")) 
        context["object_mol"]=m
    except Exception as err:
        logger.error(err)
        messages.error(req, f'{object_.pk} mol not exists or {err}')
        context["object_mol"]=''

    return render(req, "ddrug/drug/drug_detail_structure.html", context)

@login_required   
def iframe_url(req):
    context={}
    return render(req, "utils/ketcher/index.html")

@login_required   
def ketcher_test(req):
    context={}
    n=Chem.MolFromSmiles("CC(C)([C@@H]1C(O)=O)S[C@H]([C@@H]2NC([C@@H](c(cc3)ccc3O)N)=O)N1C2=O")
    m=Chem.MolToMolBlock(n)
    context["object_mol"]="\\n".join(m.split("\n"))
    return render(req, "utils/ketcher_test.html", context)

# --Drug View--
## 
class DrugListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Drug  
    template_name = 'ddrug/drug/drug_list.html' 
    filterset_class=Drug_filter
    model_fields=model.HEADER_FIELDS

##
class DrugCardView(DrugListView):
    template_name = 'ddrug/drug/drug_card.html'


    def get_queryset(self):
        queryset = super().get_queryset()  
        smiles_str=self.request.GET.get("substructure") or None
        similarity_threshold_str=self.request.GET.get("similarity_threshold") or None
        if similarity_threshold_str!=str(100) and smiles_str:
            config.tanimoto_threshold = int(similarity_threshold_str)/100
            queryset=get_mfp2_neighbors(smiles_str)
        elif smiles_str:
            queryset=Drug.objects.filter(smol__hassubstruct=QMOL(Value(smiles_str)))
        self.filterset = self.filterset_class(self.request.GET, queryset=queryset)
        order=self.get_order_by()
        if order:
            return self.filterset.qs.distinct().order_by(order)
        return self.filterset.qs.distinct()

    def get_context_data(self, **kwargs):
        try:
            context = super().get_context_data(**kwargs)
        # clearIMGfolder()
            for object_ in context["object_list"]:
                filepath=os.path.join(settings.STRUCTURE_FILES_DIR, f"{object_.pk}.svg") 
                # print(filepath)
                if os.path.exists(filepath):
                    continue
                else:
                    m=object_.smol
                    try:
                        molecule_to_svg(m, object_.pk)
                    except Exception as err:
                        messages.error(self.request, f'**{object_.pk} mol may not exists**')
        except Exception as err:
            context={}
            messages.error(self.request, err)
        return context    

##
@login_required
def detailDrug(req, pk):
    context={}
    object_=get_object_or_404(Drug, pk=pk)
    form=Drug_form(instance=object_)
    context["object"]=object_
    context["form"]=form
 
    return render(req, "ddrug/drug/drug_detail.html", context)

##
class DrugCreateView(SimplecreateView):
    form_class=Drug_form
    template_name='ddrug/drug/drug_c.html'
    
##
class DrugUpdateView(SimpleupdateView):
    form_class=Drug_form
    template_name='ddrug/drug/drug_u.html'
    model=Drug

# --Vitek Card--
class VitekcardListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_Card  
    template_name = 'ddrug/vitek_card/vitekcard_list.html' 
    filterset_class=Vitekcard_filter
    model_fields=model.HEADER_FIELDS
    context_list=''
    
    def get_context_data(self,  **kwargs):

        context =super().get_context_data( **kwargs)
        context['defaultcolumns1']='expiry_date'
        context['defaultcolumns2']='card_barcode'
        context['defaultindex1']='analysis_time'
        context['defaultindex2']='proc_date'
        context['defaultvalues']='instrument'
      
      
        data=list(context["object_list"].values())
        df=pd.DataFrame(data)
        try:
            table=pd.pivot_table(df, values=["instrument"] or None, index=["proc_date", "analysis_time"],
                        columns=["expiry_date","card_barcode"], aggfunc=np.sum).to_html(classes=["table-bordered"])
            context['table']=table
        except Exception as err:
            context['table']= 'table error : '+ str(err) 
        return context
    
   
    def post(self, request, *args, **kwargs ):
        queryset=self.get_queryset()#
        if self.request.headers.get('x-requested-with') == 'XMLHttpRequest' and self.request.method == "POST":
            selected_data=request.POST.getlist("selected_data[]") or None
            values_str=request.POST.get("values") or None
            columns_str=request.POST.get("columns") or None
            index_str=request.POST.get("index") or None
            card_barcode=request.POST.get("card_barcode")
            aggfunc_name=request.POST.get("functions")
            if selected_data:
                querydata=queryset.filter(pk__in=selected_data)
            else:
                querydata=queryset.filter(card_barcode__contains=card_barcode)
                query_send=json.dumps(list(querydata.values()), cls=DjangoJSONEncoder)   
            values=values_str or None # pivottable values
            if values:
                try:
                    table=VITEK_Card.get_pivottable(querydata=querydata, columns_str=columns_str, index_str=index_str,aggfunc=aggfunc_name, values=values)
                    response = HttpResponse(content_type='text/csv')
                    response['Content-Disposition'] = 'attachment; filename=pivottable.csv'
                    table_html=table.head().to_html(classes=["table-bordered",])
                    table_csv=table.to_csv()
                    return JsonResponse({"table_html":table_html, "table_csv":table_csv})
                except Exception as err:
                    error_message=str(err)
                    print(err)
                    return JsonResponse({"table_html":error_message,})
        return JsonResponse({})

##
@login_required
def detailVitekcard(req, pk):
    context={}
    object_=get_object_or_404(VITEK_Card, pk=pk)
    context["object"]=object_
    context["vitekid_obj"]=VITEK_ID.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    context["vitekid_fields"]=VITEK_ID.get_fields(fields=VITEK_ID.HEADER_FIELDS)
    context["vitekast_obj"]=VITEK_AST.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    context["vitekast_obj_count"]=VITEK_AST.objects.filter(card_barcode=object_.pk, astatus__gte=0).count()
    context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEK_AST.HEADER_FIELDS)

    return render(req, "ddrug/vitek_card/vitekcard_detail.html", context)

# --Vitek Ast--
## Query Across Tables:
## Card, Drug, OrgBatch, Organism
class VitekastListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_AST  
    template_name = 'ddrug/vitek_ast/vitekast_list.html' 
    filterset_class=Vitekast_filter
    model_fields=model.HEADER_FIELDS
    context_list=''

    def get_queryset(self):
        # Get the queryset however you usually would.  For example:
        queryset = super().get_queryset()

        # queryset=Vitek_CARD.objects.all().values('drug_id__drug_name','drug_id__drug_codes', 'card_barcode__orgbatch_id__organism_id__organism_name',)
        queryset=queryset.values('pk', 'drug_id__drug_name','drug_id__drug_codes', 'card_barcode__orgbatch_id__organism_id__organism_name',)
        # Then use the query parameters and the queryset to
        # instantiate a filterset and save it as an attribute
        # on the view instance for later.
        self.filterset = self.filterset_class(self.request.GET, queryset=queryset)
        # Return the filtered queryset
        order=self.get_order_by()
        if order:
            return self.filterset.qs.distinct().order_by(order)
        return self.filterset.qs.distinct()

      

## -----------
class MIC_COADDListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=MIC_COADD  
    template_name = 'ddrug/mic_coadd/mic_coadd_list.html' 
    filterset_class=MIC_COADDfilter
    model_fields=model.HEADER_FIELDS

## -------------
class MIC_COADDCardView(MIC_COADDListView):
    template_name = 'ddrug/mic_coadd/mic_coadd_card.html'
  

## -----------
class MIC_PubListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=MIC_Pub  
    template_name = 'ddrug/mic_pub/mic_pub_list.html' 
    filterset_class=MIC_Pubfilter
    model_fields=model.HEADER_FIELDS

## -------------
class MIC_PubCardView(MIC_PubListView):
    template_name = 'ddrug/mic_pub/mic_pub_card.html'

# --API Views--
## Drug
class API_Drug_List(API_FilteredListView):
    queryset = Drug.objects.all()
    serializer_class = Drug_Serializer
    filterset_class= Drug_filter

## VITEK AST
class API_VITEK_ASTList(API_FilteredListView):
    queryset = VITEK_AST.objects.all()
    serializer_class = VITEK_ASTSerializer
    filterset_class= Vitekast_filter
 
# class VITEK_ASTCreate(generics.CreateAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer

# class VITEK_ASTUpdate(generics.RetrieveUpdateAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer

# class VITEK_ASTDelete(generics.DestroyAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer
#


class Importhandler_VITEK(Importhandler):
    pass


# --upload file view--
import threading
from asgiref.sync import async_to_sync, sync_to_async
from django import forms
from apputil.utils.form_wizard_tools import ImportHandler_WizardView, UploadFileForm, StepForm_1, FinalizeForm
from django.shortcuts import render
from formtools.wizard.views import SessionWizardView
from django.core.files.storage import FileSystemStorage
from apputil.utils.validation_log import * 
from django.core.cache import cache
from django.views import View
# customized Form
class VitekValidation_StepForm(StepForm_1):
    orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',queryset=Organism_Batch.objects.filter(astatus__gte=0), widget=forms.Select(attrs={'class': 'form-control w-50'}), required=False, help_text='**Optional to choose a Organism_Batch ID',)
    field_order = ['orgbatch_id', 'confirm']

# for progress bar
def get_session_key(request):
    return request.session.session_key
# 
async def get_upload_progress(request):
    # session_key = f'upload_progress_{request.user}'
    session_key = await sync_to_async(get_session_key)(request)
    progress = await cache.aget(session_key) or {'processed': 0, 'file_name':"",'total': 0}
    print("show progress")
    print(progress)
    return JsonResponse(progress)
# 
# For get result
def get_vlog_session_key(request):   
    return f'valLog_{request.user}'
# 
async def fetchResult(request):
    cache_key = await sync_to_async(get_vlog_session_key)(request)
    valLog = await cache.aget(cache_key) or None
    if valLog is not None:
        return JsonResponse({
            'results_ready': True,
            'validation_result': valLog,
        })
    else:
        return JsonResponse({'results_ready': False})
#            
class Import_VitekView(ImportHandler_WizardView):
    
    name_step1="Validation" # step label in template
    # define more steps name
    #... 
    # define each step's form
    form_list = [
        ('upload_file', UploadFileForm),
        ('step1', VitekValidation_StepForm),
        # add more step -> StepForm
        ('finalize', FinalizeForm),
    ]
    # define template
    template_name = 'ddrug/importhandler_vitek.html'
    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filelist=[]
        self.orgbatch_id=None
        
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request
        session_key=f'upload_progress_{request.user}'

        if current_step == 'upload_file':
            print("upload and parse")
            DirName = file_location(instance=request.user)  # define file store path during file process
            files = []
            result_table=[]
            if form.is_valid():
                if 'upload_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('upload_file-multi_files'))          
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=DirName)
                    filename = fs.save(f.name, f)
                    self.filelist.append(filename)

                # Parse PDF 
                thread = threading.Thread(target=upload_VitekPDF_List, args=(request, session_key, DirName, self.filelist ), kwargs={'storage':self.storage, 'OrgBatchID': self.orgbatch_id, 'upload': False})
                thread.start()
              
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['DirName'] = DirName          
            else:
                return render(request, 'ddrug/importhandler_vitek.html', context)

        elif current_step == 'step1': # first validation
            request.session[session_key] = {'processed': 0, 'total': 0}
            upload=False#self.storage.extra_data['Confirm_to_Save']
            print("step validation again")
            form = VitekValidation_StepForm(request.POST)
            self.organism_batch=request.POST.get("upload_file-orgbatch_id") #get organism_batch
            result_table=[] # first validation result tables
            DirName=self.storage.extra_data['DirName'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name   
            thread = threading.Thread(target=upload_VitekPDF_List, args=(request, session_key, DirName, self.filelist ), kwargs={'storage':self.storage, 'OrgBatchID': self.orgbatch_id, 'upload': False})
            thread.start()         
                  
        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        print("Finalize")
        # Redirect to the desired page after finishing
        # delete uploaded files

        cache_key = f'valLog_{self.request.user}'

        filelist=self.storage.extra_data['filelist']
        for f in filelist:
            self.delete_file(f)
            print(filelist)
        cache.delete(cache_key)
        cache.delete(self.request.session.session_key)
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        # save information to context,
        # then display in templates  
        context = super().get_context_data(form=form, **kwargs)
        context['step1']=self.name_step1
        current_step = self.steps.current        
        return context

    