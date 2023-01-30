import json
import os
import asyncio
from rdkit import Chem
from django_filters.views import FilterView
import pandas as pd
import numpy as np
from django.core.serializers.json import DjangoJSONEncoder

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
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView, TemplateView
from django.utils.functional import SimpleLazyObject
from asgiref.sync import sync_to_async

from apputil.models import Dictionary, ApplicationUser
from apputil.utils import FilteredListView, get_filewithpath
from apputil.utils_dataimport import Importhandler
from apputil.views import permission_not_granted
from adjcoadd.constants import *
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from .utils import Drug_filter, Vitekcard_filter, molecule_to_svg, clearIMGfolder
from .forms import Drug_form
from .Vitek import *
   
# async_function = sync_to_async(Importhandler_VITEK.delete_file, thread_sensitive=False)
# async_function = sync_to_async(Importhandler_VITEK.post, thread_sensitive=False)
# async_function = sync_to_async(VitekcardListView.get_context_data, thread_sensitive=False)
# async_function = sync_to_async(VitekcardListView.post, thread_sensitive=False)
          
# #############################Drug View############################################
# ==========List View================================Read===========================================
class DrugListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Drug  
    template_name = 'ddrug/drug/drug_list.html' 
    filterset_class=Drug_filter
    model_fields=DRUG_FIELDs

# =============================Card View=====================================
    # editable graphic , molblock, 3D, py3Dmol 
 
class DrugCardView(DrugListView):
    template_name = 'ddrug/drug/drug_card.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # clearIMGfolder()
        for object_ in context["object_list"]:
            filepath=get_filewithpath(file_name=object_.pk)
            if os.path.exists(filepath):
                return context
            else:
                print("generate img")
                m=Chem.MolFromSmiles(object_.smiles)
                molecule_to_svg(m, object_.pk)
        return context
# ===========Detail View=============================Read============================================

# ====================================================Create===========================================
# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createDrug(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Drug_form
    if req.method=='POST':
        form=Drug_form(req.POST)
        if form.is_valid():
            instance=form.save(commit=False)
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'ddrug/drug/drug_c.html', {'form':form})
    
# ====================================================Update in Form===========================================
@login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateDrug(req, pk):
    object_=get_object_or_404(Drug, pk=pk)
    kwargs={}
    kwargs['user']=req.user 
    form=Drug_form(instance=object_)
    if req.method=='POST':
        form=Drug_form(req.POST, instance=object_)
        if form.is_valid():
            instance=form.save(commit=False)        
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'ddrug/drug/drug_u.html', {'form':form, 'object':object_})

# ================Vitek Card===========================================#
from apputil.utils import instance_dict, Validation_Log
from apputil.utils_dataimport import import_excel
from django.core.files.storage import FileSystemStorage
from pathlib import Path  
from django.core import serializers
def get_file(filename):
    from django.conf import settings
    filepath = os.path.join(settings.MEDIA_ROOT, 'table.csv')
    return filepath
import sys
class VitekcardListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_Card  
    template_name = 'ddrug/vitek_card/vitekcard_list.html' 
    filterset_class=Vitekcard_filter
    model_fields=VITEKCARD_FIELDs
    context_list=''
    
    # @classonlymethod
    # def as_view(cls, **initkwargs):
    #     view = super().as_view(**initkwargs)
    #     view._is_coroutine = asyncio.coroutines._is_coroutine
    #     return view

    # @sync_to_async
    # def get(self, request, *args, **kwargs):
      
    #     return super().get(request, *args, **kwargs)

    # @sync_to_async
    def get_context_data(self,  **kwargs):
        # get data:

        context =super().get_context_data( **kwargs)
       
        context['defaultcolumns1']='expiry_date'
        context['defaultcolumns2']='card_barcode'
        context['defaultindex1']='analysis_time'
        context['defaultindex2']='proc_date'
        context['defaultvalues']='instrument'
      
        data=list(context["object_list"].values())
        df=pd.DataFrame(data)
        table=pd.pivot_table(df, values=["instrument"], index=["proc_date", "analysis_time"],
                        columns=["expiry_date","card_barcode"], aggfunc=np.sum).to_html(classes=["table-bordered"])
        context['table']=table
    
        return context
    
    
    # @sync_to_async
    def post(self, request, *args, **kwargs ):
        queryset=self.get_queryset()#
        # receive data
        if self.request.headers.get('x-requested-with') == 'XMLHttpRequest' and self.request.method == "POST":
            selected_data=request.POST.getlist("selected_data[]") or None
            values_str=request.POST.get("values") or None
            columns_str=request.POST.get("columns") or None
            index_str=request.POST.get("index") or None
            card_barcode=request.POST.get("card_barcode")
            aggfunc_name=request.POST.get("functions")
     
            if selected_data:
                querydata=queryset.filter(pk__in=selected_data)
                print(querydata.values())
            else:
                querydata=queryset.filter(card_barcode__contains=card_barcode)
                query_send=json.dumps(list(querydata.values()), cls=DjangoJSONEncoder) 
            # get table values  
            values=values_str or None
            if values:
                try:
                    table=VITEK_Card.get_pivottable(querydata=querydata, columns_str=columns_str, index_str=index_str,aggfunc=aggfunc_name, values=values)
                    if querydata.count() >1000:
                        response = HttpResponse(content_type='text/csv')
                        response['Content-Disposition'] = 'attachment; filename=pivottable.csv'
                        table_csv=table.to_csv()
                        return JsonResponse({"table":table_csv, "msg":"large datatable automatically convert to .csv file to download...",  "table_tofront":query_send},)
                    else:
                        table_html=table.to_html(classes=["table-bordered",])
                        return JsonResponse({"table":table_html, "msg":None})
                except Exception as err:
                    error_message=str(err)
                    return JsonResponse({"table":error_message,})
        return JsonResponse({})


# ==============Vitek Card Detail===================================#
@login_required
def detailVitekcard(req, pk):
    context={}
    object_=get_object_or_404(VITEK_Card, pk=pk)
    context["object"]=object_
    context["vitekid_obj"]=VITEK_ID.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    context["vitekid_fields"]=VITEK_ID.get_fields(fields=VITEKID_FIELDs)
    context["vitekast_obj"]=VITEK_AST.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEKAST_FIELDs)

    return render(req, "ddrug/vitek_card/vitekcard_detail.html", context)


# ---------------------upload file view----------------------------
from django.conf import settings
class Importhandler_VITEK(Importhandler):

    success_url="/import-VITEK/"
    template_name='ddrug/importhandler_vitek.html'
    upload_model_type="Vitek"
    lCards={}
    lID={}
    lAst={}
    
    
    def post(self, request):
        
        dirname=settings.MEDIA_ROOT
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form
        vCard=[]
        vID=[]
        vAst=[]
        kwargs={}
        kwargs['user']=request.user
        
        self.data_model=request.POST.get('file_data')
        self.log_process=self.data_model
        myfiles=request.FILES.getlist('file_field')
        self.file_url=[]
        try:
            if form.is_valid():
                self.lCards.clear()
                self.lID.clear()
                self.lAst.clear()
                self.table_name.clear()
                self.validate_result.clear()
                self.file_report.clear()
                print(myfiles)
                print(self.file_report)
                # scan_results = cd.instream(myfile) # scan_results['stream'][0] == 'OK' or 'FOUND'
                for f in myfiles:
                    fs=FileSystemStorage()
                    print(f"fs is {fs}")
                    filename=fs.save(str(request.user)+"_"+f.name, f)
                    
                    try:
                        vCards, vID, vAst=process_VitekPDF(DirName=dirname, PdfName=filename)#import_excel(fs.url(filename), self.data_model)
                        print("file checked")
                        self.lCards[filename]=vCards
                        self.lID[filename]=vID
                        self.lAst[filename]=vAst
                        self.file_url.append(fs.url(filename))
                    except Exception as err:
                        self.delete_file(file_path=fs.url(filename))
                        messages.warning(request, f'{filename} contains {err} erro , cannot to upload, choose a correct file')
                        
                   
                context['file_pathlist']=self.file_url
                context['data_model']=self.data_model
                context['cards']=self.lCards
                context['ids']=self.lID    

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')

        if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST": 
         

            process_name=request.POST.get('type')
            file_pathlist=request.POST.getlist("filepathlist[]")
            print(f'selected : {file_pathlist}')
            # for f in file_pathlist:
            print(f"cardlist is: {self.lCards}")
            data_model=request.POST.get("datamodel")
          
            if process_name=='Validation':
                if any(self.lCards.values()):
                    for key in self.lCards:
                        self.table_name.append(key)
                      
                        for e in self.lCards[key]:
                            djCard=VITEK_Card.check_from_dict(e, self.vLog)
                            
                            if self.validate_result.get(key, None):
                                self.validate_result[key].append("Card-"+str(djCard.validStatus))
                            else:
                                self.validate_result[key]=[("Card-"+str(djCard.validStatus))]
                            
                                # result_report["key"]=["Card-"+str(djCard.validStatus)]
                            
                            self.file_report.append(key+str(self.vLog.show()))
                
                if any(self.lID.values()):
                    for key in self.lID:
                        self.table_name.append(key)
                        for e in self.lID[key]:
                            djID=VITEK_ID.check_from_dict(e, self.vLog)
                            if self.validate_result.get(key, None):
                                self.validate_result[key].append("ID-"+str(djID.validStatus))
                            else:
                                self.validate_result[key]=[("ID-"+str(djID.validStatus))]
                            # self.validate_result.append(f"{key}+ID status: {djID.validStatus}")
                            self.file_report.append(str(self.vLog.show()))
                
                if any(self.lAst.values()):
                    for key in self.lAst:
                        self.table_name.append(key)
                        for e in self.lAst[key]:
                            djAst=VITEK_AST.check_from_dict(e, self.vLog)
                            if self.validate_result.get(key, None):
                                self.validate_result[key].append("Ast-"+str(djAst.validStatus))
                            else:
                                self.validate_result[key]=[("Ast-"+str(djAst.validStatus))]
                            # self.validate_result.append(f"{key}+Ast status: {djAst.validStatus}")
                            self.file_report.append(key+str(self.vLog.show()))

                       
        
                return JsonResponse({"table_name":" ||Vitek| ".join(self.table_name), 'validate_result':str(self.validate_result), 'file_report':"," .join(self.table_name)})                                   
       
            elif process_name=='Cancel':
                # Cancel Task
                self.file_url=[]
                for f in file_pathlist:
                    self.delete_file(file_path=f)
                    self.file_url.append(f)
                self.lCards.clear()
                self.lID.clear()
                self.lAst.clear()
                self.table_name.clear()
                self.validate_result.clear()
                self.file_report.clear()
                return JsonResponse({"table_name":(",").join( self.table_name), 'validate_result':"Deleted", 'file_report':f"delete files: {self.file_url}"})                                   
    
            elif process_name=='DB_Validation':
                print("start saving to db")
                print(self.lCards)

                if any(self.lCards.values()):
                    for key in self.lCards:
                        self.table_name.append(key)
                        for e in self.lCards[key]:
                            djCard=VITEK_Card.check_from_dict(e, self.vLog)
                            if djCard.validStatus:
                                try:
                                    djCard.save(**kwargs)
                                except Exception as err:
                                    self.validate_result[key].append(f"catch Exception Card {err}")
                            self.validate_result[key].append(f"CARD status: {djCard.validStatus}")
                            self.file_report.append(str(self.vLog.show()))
                
                if any(self.lID.values()):
                    for key in self.lID:
                        self.table_name.append(key)
                        for e in self.lID[key]:
                            djID=VITEK_ID.check_from_dict(e, self.vLog)
                            if djID.validStatus:
                                try:
                                    djID.save(**kwargs)
                                except Exception as err:
                                    self.validate_result[key].append(f"catch Exception ID {err}")
                            self.validate_result[key].append(f"ID status: {djID.validStatus}")
                            self.file_report.append(str(self.vLog.show()))
                
                if any(self.lAst.values()):
                    for key in self.lAst:
                        self.table_name.append(key)
                        for e in self.lAst[key]:
                            djAst=VITEK_AST.check_from_dict(e, self.vLog)
                            if djAst.validStatus:
                                try:
                                    djAst.save(**kwargs)
                                except Exception as err:
                                    self.validate_result[key].append(f"catch Exception Ast {err}")
                            self.validate_result[key].append(f"Ast status: {djAst.validStatus}")
                            self.file_report.append(str(self.vLog.show()))
                
                print(f'lcards clear is {self.lCards}')
                for f in file_pathlist:
                    self.delete_file(file_path=f)
                               
                return JsonResponse({"table_name":" ||Vitek| ".join(self.table_name), 'validate_result':str(self.validate_result), 'file_report':(",").join(self.file_report), 'status':"Data Saved! uploaded files clear!"})

           
        return render(request, 'ddrug/importhandler_vitek.html', context)




      


  
