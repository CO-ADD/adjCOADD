import os
import json
from rdkit import Chem
from django_filters.views import FilterView

from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction, IntegrityError
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.utils.functional import SimpleLazyObject
from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import FilteredListView
from apputil.utils.views_base import permission_not_granted, SimplecreateView, SimpleupdateView,  SimpledeleteView
from adjcoadd.constants import *
from .models import  Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from .forms import (CreateOrganism_form, UpdateOrganism_form, Taxonomy_form, 
                    Batch_form, Batchupdate_form, Stock_createform, Stock_form, Culture_form, Cultureupdate_form,
                    Organismfilter, Taxonomyfilter, Batchfilter, Stockfilter)

from ddrug.models import VITEK_AST, MIC_COADD
from .utils.data_visual import data_frame_style, pivottable_style

          
# --TAXONOMY Views--
##
class TaxonomyListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Taxonomy  
    template_name = 'dorganism/taxonomy/taxonomy_list.html' 
    filterset_class=Taxonomyfilter
    model_fields=model.HEADER_FIELDS

##
class TaxonomyCardView(TaxonomyListView):
    template_name = 'dorganism/taxonomy/taxonomy_card.html'
    
##
@login_required
def detailTaxonomy(req, slug=None):
    context={}
    object_=get_object_or_404(Taxonomy, urlname=slug)
    context["object"]=object_
    context['form']=Taxonomy_form(instance=object_)
    return render(req, "dorganism/taxonomy/taxonomy_detail.html", context)

##
class TaxonomyCreateView(SimplecreateView):
    form_class=Taxonomy_form
    template_name='dorganism/taxonomy/taxonomy_c.html'
    transaction_use = 'dorganism'
        
##
class TaxonomyUpdateView(SimpleupdateView):
    form_class=Taxonomy_form
    template_name='dorganism/taxonomy/taxonomy_u.html'
    model=Taxonomy
    transaction_use = 'dorganism'

##
class TaxonomyDeleteView(SimpledeleteView):
    model = Taxonomy
    transaction_use = 'dorganism'


#--ORGANISM Views--
##
class OrganismListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = Organism  
    template_name = 'dorganism/organism/organism_list.html'
    filterset_class = Organismfilter
    model_fields = model.HEADER_FIELDS
    
##  
class OrganismCardView(OrganismListView):
    template_name = 'dorganism/organism/organism_card.html'
    model = Organism  
    model_fields = model.CARDS_FIELDS
    

##
    # =============================step 2. Create new record by form===================#
@login_required
def createOrganism(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    form=CreateOrganism_form()
    if req.method=='POST':
        Organism_Name=req.POST.get('search_organism') # -1. Ajax Call(/search_organism/) find Foreignkey orgname
        form=CreateOrganism_form( Organism_Name, req.POST,) # -2. get create-form with ajax call result
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'): # -3. write new entry in atomic transaction
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, f'create failed due to {form.errors} error')
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dorganism/organism/organism_c.html', { 'form':form, }) 


##
@login_required
def detailOrganism(request, pk):
    """
    - Detail view handle organism single entry
    display,update and delete.
    - related table overview display, update, create and delete.
    - related table are: batch, stock, culture.
    - data visual table: dataframe and pivot- table
    """
   
    from django.db.models import Count
    from apputil.forms import Image_form, Document_form
    context={}
    object_=get_object_or_404(Organism, organism_id=pk)
    try:
        form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel,}, instance=object_)
    except Exception as err:
        print(err)
    context["object"]=object_
    context["form"]=form
    context["image_form"]=Image_form
    context["doc_form"]=Document_form

    # data in related tables
    context["batch_obj"]=Organism_Batch.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["batch_obj_count"]=context["batch_obj"].count() if context["batch_obj"].count()!=0 else None
    context["batch_fields"]=Organism_Batch.get_fields()
    context["stock_count"]=Organism_Batch.objects.annotate(number_of_stocks=Count('orgbatch_id')) 
    context["cultr_obj"]=Organism_Culture.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["cultr_obj_count"]=context["cultr_obj"].count() if context["cultr_obj"].count()!=0 else None
    context["cultr_fields"]=Organism_Culture.get_fields() 
    if 'organism_id' in context["cultr_fields"]:
        context["cultr_fields"].remove('organism_id')    # customize HEADER_FIELDS
    context["vitekast_obj"]=SimpleLazyObject(lambda: VITEK_AST.objects.filter(organism=object_.organism_name, astatus__gte=0))
    context["vitekast_obj_count"]=context["vitekast_obj"].count() if context["vitekast_obj"].count()!=0 else None
    context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEK_AST.HEADER_FIELDS)

    # data in pivotted and highlighted Tables
    if request.method == 'POST':
        displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source']
        context["table"] = data_frame_style(pk, displaycols)['style_table']
        context["df_entries"] = data_frame_style(pk, displaycols)['df_entries']
        context["pivottable"] = pivottable_style(pk)
        return render(request, "dorganism/organism/organism_mic.html", context)
    


    return render(request, "dorganism/organism/organism_detail.html", context)

##
@login_required
def updateOrganism(req, pk):
    object_=get_object_or_404(Organism, organism_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel, 'assoc_images': [i.image_file for i in object_.assoc_images.all()], 'assoc_documents': [i.doc_file for i in object_.assoc_documents.all()]}, instance=object_)
    if object_.organism_name.org_class: # Organism_Class_str for display class
        Organism_Class_str=object_.organism_name.org_class.dict_value
    else:
        Organism_Class_str="No Class"
    if req.method=='POST':
        try:
            with transaction.atomic(using='dorganism'):        # testing!
                obj = Organism.objects.select_for_update().get(organism_id=pk)
                #If update Organism Name
                if  req.POST.get('search_organism'):
                    Organism_Name_str=req.POST.get('search_organism')
                    Organism_new_obj=get_object_or_404(Taxonomy, organism_name=Organism_Name_str)
                    form=UpdateOrganism_form(Organism_Name_str, req.POST, instance=obj)
                    #-Not allow to update name in different class--
                    if Organism_new_obj.org_class.dict_value and Organism_new_obj.org_class.dict_value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                else:
                    form=UpdateOrganism_form(object_.organism_name, req.POST, instance=obj) 
                
                if form.is_valid():       
                    instance=form.save(commit=False)
                    instance.save(**kwargs)
                    # form.save_m2m() 
                    return redirect(req.META['HTTP_REFERER'])
                else:
                    messages.warning(req, f'Update failed due to {form.errors} error')
                   
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])
  
    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dorganism/organism/organism_u.html", context)

##

class OrganismDeleteView(SimpledeleteView):
    model = Organism
    transaction_use = 'dorganism'


   
#--Batch Views--
##
class BatchCardView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism_Batch 
    template_name = 'dorganism/organism/batch/batch_card.html' 
    filterset_class=Batchfilter
    model_fields=model.HEADER_FIELDS

## 
@login_required
def createBatch(req, organism_id):
    kwargs={'user': req.user}
    form=Batch_form()
    
    if req.method=='POST':
        form=Batch_form(req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.organism_id=get_object_or_404(Organism, pk=organism_id)              
                    instance.save(**kwargs)
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/organism/batch/batch_c.html', { 'form':form, 'organism_id':organism_id}) 

## here used HTMX
from apputil.utils.views_base import HtmxupdateView
class BatchUpdateView(HtmxupdateView):
    form_class=Batchupdate_form
    template_name="dorganism/organism/batch/batch_u.html"
    template_partial="dorganism/organism/batch/batch_tr.html"
    model=Organism_Batch
    transaction_use = 'dorganism'

##
class BatchDeleteView(SimpledeleteView):
    model = Organism_Batch
    transaction_use = 'dorganism'

# --Stock Views--
# view in organism detail views
## here is response to an Ajax call
## to send data to child datatable 
@user_passes_test(lambda u: u.has_permission('Read'), login_url='permission_not_granted') 
def stockList(req, pk):
    res=None
    if req.method == 'GET':
        batch_id=req.GET.get('Batch_id')
        object_=get_object_or_404(Organism_Batch, orgbatch_id=batch_id)#Organism_Batch.objects.get(orgbatch_id=batch_id)
        qs=OrgBatch_Stock.objects.filter(orgbatch_id=object_, astatus__gte=0, n_left__gt=1) # n_left show when bigger or equal to 2
        data=[]
        for i in qs:
            item={
                "stock_id":i.pk,
                "stock_type": str(i.stock_type.dict_value) or 'no data',
                "location_freezer":str(i.location_freezer) or 'no data',
                "location_rack": str(i.location_rack),
                "location_col": str(i.location_column),
                "location_slot": str(i.location_slot),
                "stock_date": str(i.stock_date.strftime("%d-%m-%Y")),
                "n_left": str(i.n_left) if i.n_left >=1 else " ",
                "n_created": str(i.n_created),
                "stock_notes":str(i.stock_note),
                "biologist": str(i.biologist),
            }
            data.append(item)          
        res=data        
        return JsonResponse({'data':res})
    return JsonResponse({})
#
# Overview Stocks
##
# from apputil.utils.flex_pivottable import flex_pivottable 
class StockListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = OrgBatch_Stock  
    template_name = 'dorganism/organism/batch_stock/stock_list.html'
    filterset_class = Stockfilter
    model_fields = model.HEADER_FIELDS 

    def get_context_data(self,  **kwargs):

        context =super().get_context_data( **kwargs)
        print(self.request.GET)

        # context['defaultcolumns1']='expiry_date'
        # context['defaultcolumns2']='card_barcode'
        # context['defaultindex1']='analysis_time'
        # context['defaultindex2']='proc_date'
        # context['defaultvalues']='instrument'
        
        # context['pivottable_stock'] = flex_pivottable(self.request, model=self.model)
        # print(context['pivottable_stock'])
        return context

      


@login_required
def createStock(req, orgbatch_id):
    kwargs={}
    kwargs['user']=req.user
    form=Stock_createform() #Stock_createform(initial={"orgbatch_id":orgbatch_id},)
    if req.method=='POST':
        form=Stock_createform(req.POST)
        if form.is_valid():
            orgbatch_id=orgbatch_id
            stock_type=req.POST.get("stock_type")
            stock_date=req.POST.get("stock_date")
            n_created=req.POST.get("n_created")
            kwargs['orgbatch_id']=orgbatch_id
            kwargs['stock_type']=stock_type
            kwargs['stock_date']=stock_date
            kwargs['n_created']=n_created

            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'wrong {form.errors}')
    return render(req, 'dorganism/organism/batch_stock/stock_c.html', { 'form':form, 'orgbatch_id':orgbatch_id }) 

##
@login_required
def updateStock(req, pk):
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Stock_form(instance=object_)
    if req.method == 'POST' and req.headers.get('x-requested-with') == 'XMLHttpRequest':
        # process the data sent in the AJAX request
        n_left_value=req.POST.get('value')
        object_.n_left=int(n_left_value)-1
        object_.save(**kwargs)
        response_data = {'result': str(object_.n_left)}
        return JsonResponse(response_data)
    if req.method=='POST':
        form=Stock_form(req.POST, instance=object_)
        if "cancel" in req.POST:
            return redirect(req.META['HTTP_REFERER'])
        else:
            try:
                with transaction.atomic(using='dorganism'):      
                    obj = OrgBatch_Stock.objects.select_for_update().get(pk=pk)
                    try:
                        if form.is_valid():               
                            instance=form.save(commit=False)
                            instance.save(**kwargs)
                            return redirect(req.META['HTTP_REFERER'])
                            
                    except Exception as err:
                        print(f'form erroro is {form.errors} and error {err}')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
                return redirect(req.META['HTTP_REFERER'])
    context={
        "form":form,
        "object":object_,
    }
    return render(req, "dorganism/organism/batch_stock/stock_u.html", context)

##
class StockDeleteView(SimpledeleteView):
    model = OrgBatch_Stock
    transaction_use = 'dorganism'

# --Culture Views--
## 
@login_required
def createCulture(req, organism_id):
    print(f"id is {organism_id}")
    kwargs={'user' : req.user}
    form=Culture_form()
    print('form')
    if req.method=='POST':
        print("save") 
        form=Culture_form(req.POST)
        if form.is_valid():
            print("save") 
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False)
                    instance.organism_id=get_object_or_404(Organism, pk=organism_id)
                    print(instance.organism_id) 
                    instance.save(**kwargs)                  
                    # message=instance.save(**kwargs)
                    # if type(message)==Exception:
                    #     messages.error(req, f'{message} happes')
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(form.errors)
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/organism/culture/culture_c.html', { 'form':form, 'organism_id':organism_id}) 

## Here used HTMX
class CultureUpdateView(HtmxupdateView):
    form_class=Cultureupdate_form
    template_name="dorganism/organism/culture/culture_u.html"
    template_partial="dorganism/organism/culture/culture_tr.html"
    model=Organism_Culture
    transaction_use = 'dorganism'

##
class CultureDeleteView(SimpledeleteView):
    model = Organism_Culture
    transaction_use = 'dorganism'

