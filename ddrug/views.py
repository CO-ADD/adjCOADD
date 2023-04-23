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
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub
from .utils import (molecule_to_svg, 
                    clearIMGfolder, get_mfp2_neighbors)
from .forms import Drug_form, Drug_filter, Vitekcard_filter, Vitekast_filter, MIC_COADDfilter, MIC_Pubfilter
from .serializers import Drug_Serializer, VITEK_ASTSerializer
from .util_vitek import *

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
    print(context["object_mol"])
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
                        print(object_.pk)
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
        print(context['fields'])
      
        data=list(context["object_list"].values())
        print(data)
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
            print(f'columns_str is {columns_str}')
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
        print(f'query: {queryset}')
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


# --upload file view--

class Importhandler_VITEK(Importhandler):

    success_url="/import-VITEK/"
    template_name='ddrug/importhandler_vitek.html'
    upload_model_type="Vitek"
    lCards={}   # self.lCards, lID and lAst store results parsed by all uploaded files with key-filename, value-parsed result array
    lID={}
    lAst={}
    # vLog = Validation_Log("Vitek-pdf")
    
    def post(self, request, process_name):
        location=file_location(request) # define file store path during file process
        form = self.form_class(request.POST, request.FILES)
        context, kwargs = {}, {}
        context['form'] = form
        kwargs['user']=request.user
        vCard, vID, vAst=[], [], []   # process variables for store parse result from single file
        
        self.file_url=[]
        self.data_model=request.POST.get('file_data') or None
        myfiles=request.FILES.getlist('file_field')
      
        try:
        # Uploading Verifying     
            if form.is_valid():
                self.lCards.clear()
                self.lID.clear()
                self.lAst.clear()
                
        # Uploading Parsing 
                for f in myfiles:
                    fs=OverwriteStorage(location=location)
                    filename=fs.save(f.name, f)
                    print(f"fs is {fs} location is {location}, filename is {filename}")
              
                    try:
                        vCards, vID, vAst=process_VitekPDF(DirName=location, PdfName=filename)
                        print("file checked")
                        self.lCards[filename]=vCards
                        self.lID[filename]=vID
                        self.lAst[filename]=vAst
                        self.file_url.append(filename)
                    except Exception as err:
                        self.delete_file(file_name=filename)
                        messages.warning(request, f'{filename} contains {err} erro , cannot to upload, choose a correct file')
                        return render(request, 'ddrug/importhandler_vitek.html', context)
                  
                context['file_list']=self.file_url
                context['data_model']=self.data_model
                context['cards']=self.lCards
                context['ids']=self.lID   
                print(self.lAst) 

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')

        # Ajax Call steps for validation, cancel/delete or savetoDB 
        if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":          
            vLog = Validation_Log("Vitek-pdf") # create process Log instance
            process_name=request.POST.get('type') # get process type: Validation, Delete, DB_Validation(Save to DB)
            file_list=request.POST.getlist("select_file_list[]") #get selected fileName List
            data_model=request.POST.get("datamodel")
            self.validate_result.clear()
            self.file_report.clear()
            lCards, lID, lAst={},{},{}
            for f in file_list:
                lCards[f], lID[f], lAst[f]=self.lCards[f], self.lID[f], self.lAst[f]     
        # Validating
            if process_name=='Validation':        
                self.validates(lCards, VITEK_Card, vLog, self.validate_result, self.file_report, save=False, **kwargs)
                self.validates(lID, VITEK_ID, vLog, self.validate_result, self.file_report, save=False, **kwargs)
                self.validates(lAst, VITEK_AST, vLog, self.validate_result, self.file_report, save=False, **kwargs)
               
                return JsonResponse({ 'validate_result':str(self.validate_result), 'file_report':str(self.file_report).replace("\\", "").replace("_[", "_").replace("]_", "_"), 'status':'validating'})                                   
        # Cancel and deleting Task                
            elif process_name=='Delete':
                for f in file_list:
                    try:
                        self.delete_file(file_name=f)                   
                    except Exception as err:
                        return JsonResponse({"status":"Delete", "systemErr":"File not existed!"})               
                return JsonResponse({"status":"Delete"})                                   
        # Saving
            elif process_name=='DB_Validation':
                print("start saving to db")                         
                self.validates(lCards, VITEK_Card, vLog, self.validate_result, self.file_report, save=True, **kwargs)
                self.validates(lID, VITEK_ID, vLog, self.validate_result, self.file_report, save=True, **kwargs)
                self.validates(lAst, VITEK_AST, vLog, self.validate_result, self.file_report, save=True, **kwargs)                     
                return JsonResponse({ 'validate_result':str(self.validate_result), 'file_report':str(self.file_report).replace("\\", "").replace("_[", "_").replace("]_", "_"), 
                'status':'SavetoDB', 'savefile':str(file_list)})
        context["process_name"]=process_name   
        return render(request, 'ddrug/importhandler_vitek.html', context)


