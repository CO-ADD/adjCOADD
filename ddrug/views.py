import json
import os
from rdkit import Chem
from django_filters.views import FilterView
import pandas as pd
import numpy as np
from django.core.serializers.json import DjangoJSONEncoder
from time import localtime, strftime

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


from apputil.models import Dictionary, ApplicationUser
from apputil.utils import FilteredListView, get_filewithpath, file_location
from apputil.utils_dataimport import Importhandler
from apputil.views import permission_not_granted
from adjcoadd.constants import *
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from .utils import Drug_filter, Vitekcard_filter, molecule_to_svg, clearIMGfolder
from .forms import Drug_form
from .Vitek import *
   

          
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
from apputil.utils import instance_dict, Validation_Log, OverwriteStorage
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
            print(f'columns_str is {columns_str}')
            if selected_data:
                querydata=queryset.filter(pk__in=selected_data)
            else:
                querydata=queryset.filter(card_barcode__contains=card_barcode)
                query_send=json.dumps(list(querydata.values()), cls=DjangoJSONEncoder) 
            # get table values  
            values=values_str or None
            if values:
                try:
                    table=VITEK_Card.get_pivottable(querydata=querydata, columns_str=columns_str, index_str=index_str,aggfunc=aggfunc_name, values=values)
                    # if querydata.count() >1000:
                    response = HttpResponse(content_type='text/csv')
                    response['Content-Disposition'] = 'attachment; filename=pivottable.csv'
                    table_html=table.head().to_html(classes=["table-bordered",])
                    table_csv=table.to_csv()
                    return JsonResponse({"table_html":table_html, "table_csv":table_csv, "table_tofront":query_send},)
                    # else:
                    #     table_json=table.to_json()
                    #     return JsonResponse({"table":table_html, "msg":None, "table_json":table_json})
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
    # vLog = Validation_Log("Vitek-pdf")
    
    def post(self, request):
                
        location=file_location(request)
        form = self.form_class(request.POST, request.FILES)
        context = {}
        context['form'] = form
        vCard=[]
        vID=[]
        vAst=[]
        kwargs={}
        kwargs['user']=request.user
        
        self.file_url=[]
        self.data_model=request.POST.get('file_data')
        # self.log_process=self.data_model
        myfiles=request.FILES.getlist('file_field')
      
        try:
            if form.is_valid():
                self.lCards.clear()
                self.lID.clear()
                self.lAst.clear()
                
                # scan_results = cd.instream(myfile) # scan_results['stream'][0] == 'OK' or 'FOUND'
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
                        
                   
                context['file_pathlist']=self.file_url
                context['data_model']=self.data_model
                context['cards']=self.lCards
                context['ids']=self.lID    

        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')

        if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":          
            # steps for validation, cancel/delete or savetoDB 
            vLog = Validation_Log("Vitek-pdf")
            print(self.file_url)
            process_name=request.POST.get('type')
            file_pathlist=request.POST.getlist("filepathlist[]")
            data_model=request.POST.get("datamodel")
            # self.validatefile_name.clear()
            self.validate_result.clear()
            self.file_report.clear()
            self.lCards.clear()
            self.lID.clear()
            self.lAst.clear()
        
            if process_name=='Validation':
                for f in file_pathlist:
                    # self.validatefile_name.append(f)
                    if self.lCards.get(f, None):
                        print(self.lCards)
                        pass
                    else:
                        try:
                            vCards, vID, vAst=process_VitekPDF(DirName=location, PdfName=f)
                            print("file checked")
                            self.lCards[f]=vCards
                            self.lID[f]=vID
                            self.lAst[f]=vAst    
                        except Exception as err:
                            self.delete_file(file_name=f)
                            messages.warning(request, f'{filename} contains {err} erro , cannot to upload, choose a correct file')
     
                self.validates(self.lCards, VITEK_Card, vLog, self.validate_result, self.file_report, save=False, **kwargs)
                self.validates(self.lID, VITEK_ID, vLog, self.validate_result, self.file_report, save=False, **kwargs)
                self.validates(self.lAst, VITEK_AST, vLog, self.validate_result, self.file_report, save=False, **kwargs)
       
               
                return JsonResponse({ 'validate_result':str(self.validate_result), 'file_report':str(self.file_report).replace("\\", "").replace("_[", "_").replace("]_", "_"), 'status':'validating'})                                   
       
            elif process_name=='Cancel':
                # Cancel Task
                self.file_url=[]
                for f in file_pathlist:
                    self.delete_file(file_name=f)
                    self.file_url.append(f)
                self.lCards.clear()
                self.lID.clear()
                self.lAst.clear()
                # "validatefile_name":(",").join( self.validatefile_name), 
                return JsonResponse({"status":"Delete"})                                   
    
            elif process_name=='DB_Validation':
                
                print("start saving to db")
                print(self.lCards)
                self.validates(self.lCards, VITEK_Card, vLog, self.validate_result, self.file_report, save=True, **kwargs)
                self.validates(self.lID, VITEK_ID, vLog, self.validate_result, self.file_report, save=True, **kwargs)
                self.validates(self.lAst, VITEK_AST, vLog, self.validate_result, self.file_report, save=True, **kwargs)
               
                #   "validatefile_name":" ||Vitek| ".join(self.validatefile_name),          
                return JsonResponse({ 'validate_result':str(self.validate_result), 'file_report':str(self.file_report).replace("\\", "").replace("_[", "_").replace("]_", "_"), 'status':"Data Saved! uploaded files clear!"})

           
        return render(request, 'ddrug/importhandler_vitek.html', context)




      


  
