import os
from rdkit import Chem
from django_filters.views import FilterView
import pandas as pd
import numpy as np

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
from apputil.utils import FilteredListView
from apputil.views import permission_not_granted
from adjcoadd.constants import *
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from .utils import get_filewithpath, Drug_filter, Vitekcard_filter, molecule_to_svg, clearIMGfolder
from .forms import Drug_form
   
          
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
class VitekcardListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_Card  
    template_name = 'ddrug/vitek_card/vitekcard_list.html' 
    filterset_class=Vitekcard_filter
    model_fields=VITEKCARD_FIELDs

    # def get_context_data(self, **kwargs):
    #     context = super().get_context_data(**kwargs)
    #     print(context["page_obj"])
    #     # df=pd.DataFrame()


class VitekcardPivotView(VitekcardListView):
    template_name='ddrug/vitek_card/vitekcard_pivotable.html'


    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        data=list(context["object_list"].values())
        
        values_str=self.request.GET.get("values") or None
        columns_str=self.request.GET.get("columns") or None
        indexs_str=self.request.GET.get("index") or None
        # if (indexs_str):
        values=values_str or None
        if values:
            df=pd.DataFrame(data)
            columns=columns_str.split[","] 
            index=index_str.split[","] 

            table=pd.pivot_table(df, values=values, index=index,
                    columns=columns, aggfunc=np.sum).to_html(classes=["table-bordered", "table-striped", "table-hover"]) 
        
        else:
            #example default pivottable
            data_example=list(VITEK_Card.objects.filter(proc_date='2021-06-01').values())
            df_example=pd.DataFrame(data_example)
            values_example='analysis_time' 
            columns_example=['card_type_id','card_barcode', ] 
            index_example=['acreated_id',  'orgbatch_id_id',] 
            table_example=pd.pivot_table(df_example, values=values_example, index=index_example,
                    columns=columns_example, aggfunc=np.sum).to_html(classes=["table-bordered", "table-striped", "table-hover"]) 
            
            context['columns']=columns_example
            context['table']=table_example
            return context

        context['columns']= columns
        context['table']=table
        return context
        # df=pd.DataFrame()

# ==============Vitek Card Detail===================================#
@login_required
def detailVitekcard(req, pk):
    context={}
    object_=get_object_or_404(VITEK_Card, pk=pk)
    context["object"]=object_
    context["vitekid_obj"]=VITEK_ID.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    context["vitekid_fields"]=VITEK_ID.get_fields(fields=VITEKID_FIELDs)
    context["vitekast_obj"]=VITEK_AST.objects.filter(card_barcode=object_.pk, astatus__gte=0)
    print(context["vitekast_obj"])
    context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEKAST_FIELDs)

    return render(req, "ddrug/vitek_card/vitekcard_detail.html", context)
