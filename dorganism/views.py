import os
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
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView, TemplateView
from django.utils.functional import SimpleLazyObject
from apputil.models import Dictionary, ApplicationUser
from apputil.utils import FilteredListView
from apputil.views import permission_not_granted
from adjcoadd.constants import *
from .models import  Organism, Taxonomy, Organism_Batch, OrgBatch_Stock, Organism_Culture
from .utils import  Organismfilter, Taxonomyfilter, Batchfilter
from .forms import (CreateOrganism_form, UpdateOrganism_form, Taxonomy_form, 
                    Batch_form, Batchupdate_form, Stock_createform, Stock_form, Culture_form, Cultureupdate_form)
   
          
# #############################TAXONOMY View############################################
# ==========List View================================Read===========================================
class TaxonomyListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Taxonomy  
    template_name = 'dorganism/taxonomy/taxonomy_list.html' 
    filterset_class=Taxonomyfilter
    model_fields=model.HEADER_FIELDS

# =============================Card View=====================================

class TaxonomyCardView(TaxonomyListView):
    template_name = 'dorganism/taxonomy/taxonomy_card.html'

    
# ===========Detail View=============================Read============================================
@login_required
def detailTaxonomy(req, slug=None):
    context={}
    object_=get_object_or_404(Taxonomy, urlname=slug)
    context["object"]=object_
    context['form']=Taxonomy_form(instance=object_)
    return render(req, "dorganism/taxonomy/taxonomy_detail.html", context)

# ====================================================Create===========================================
@login_required
def createTaxonomy(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Taxonomy_form
    if req.method=='POST':
        form=Taxonomy_form(req.POST)
        if form.is_valid():
            instance=form.save(commit=False)
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/taxonomy/taxonomy_c.html', {'form':form})
    
# ====================================================Update in Form===========================================
@login_required
def updateTaxonomy(req, slug=None):
    object_=get_object_or_404(Taxonomy, urlname=slug)
    kwargs={}
    kwargs['user']=req.user 
    form=Taxonomy_form(instance=object_)
    if req.method=='POST':
        form=Taxonomy_form(req.POST, instance=object_)
        if form.is_valid():
            instance=form.save(commit=False)        
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'dorganism/taxonomy/taxonomy_u.html', {'form':form, 'object':object_})

# ====================================================Delete===========================================
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted')
def deleteTaxonomy(req, slug=None):
    kwargs={}
    kwargs['user']=req.user 
    object_=get_object_or_404(Taxonomy, urlname=slug)
    try:
        if req.method=='POST':
            object_.delete(**kwargs)
    except Exception as err:
        print(err) 
    return redirect("taxo_card")

############################################### ORGANISM View ###########################################
# ==============================List View ================================
class OrganismListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism  
    template_name = 'dorganism/organism/organism_list.html'
    filterset_class=Organismfilter
    model_fields=model.HEADER_FIELDS
    
# =============================Card View=====================================
    
class OrganismCardView(OrganismListView):
    template_name = 'dorganism/organism/organism_card.html'

# ======================================================================CREATE==========================================#
    # ==Step1. Ajax Call(def search_organism in utils) search Taxonomy(for all models using Taxonomy as ForeignKey)=====#
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
        Organism_Name=req.POST.get('search_organism')
        form=CreateOrganism_form( Organism_Name, req.POST,)
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    print("org saved")
                    return redirect(req.META['HTTP_REFERER'])

            except IntegrityError as err:
                    print("error")
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
        

    return render(req, 'dorganism/organism/organism_c.html', { 'form':form, }) 


#=========================================Organism detail table with updating in detail table========================================================================================
@login_required
def detailOrganism(req, pk):
    context={}
    object_=get_object_or_404(Organism, organism_id=pk)
    form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel}, instance=object_)
    context["object"]=object_
    if object_.strain_type:
        context["strain_type"]=",".join(object_.strain_type)
    else:
        context["strain_type"]=" "
    if object_.strain_panel:
        context["strain_panel"]=",".join(object_.strain_panel)
    else:
        context["strain_panel"]=" "
    context["form"]=form
    context["batch_obj"]=Organism_Batch.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["batch_fields"]=Organism_Batch.get_fields()
    context["cultr_obj"]=Organism_Culture.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["cultr_fields"]=Organism_Culture.get_fields()

    return render(req, "dorganism/organism/organism_detail.html", context)

#======================================================================Update Organism=================================================================================
@login_required
def updateOrganism(req, pk):
    object_=get_object_or_404(Organism, organism_id=pk)
    kwargs={}
    kwargs['user']=req.user
    #This can be minimized when all organism have classes... ----------------
    form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel}, instance=object_)
    if object_.organism_name.org_class:
        Organism_Class_str=object_.organism_name.org_class.dict_value
    else:
        Organism_Class_str="No Class"
    #-------------------------------------------------------------------------
    if req.method=='POST':

        try:
            with transaction.atomic(using='dorganism'):        # testing!
                obj = Organism.objects.select_for_update().get(organism_id=pk)
                #------------------------If update Organism Name-----------------------------------
                if  req.POST.get('search_organism'):
                    Organism_Name_str=req.POST.get('search_organism')
                    Organism_new_obj=get_object_or_404(Taxonomy, organism_name=Organism_Name_str)
                    form=UpdateOrganism_form(Organism_Name_str, req.POST, instance=obj)
                    #-----------------------Not allow to update name in different class--------
                    if Organism_new_obj.org_class.dict_value and Organism_new_obj.org_class.dict_value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                    #-----------------------Not allow to update name in different class--------
                #------------------------If update Organism Name-----------------------------------
                else:
                    form=UpdateOrganism_form(object_.organism_name, req.POST, instance=obj) 
                try:
                    if form.is_valid():                  
                        instance=form.save(commit=False)
                        instance.save(**kwargs)
                        return redirect(req.META['HTTP_REFERER'])
                except Exception as err:
                    print(err)
                   
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])
  
    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dorganism/organism/organism_u.html", context)

# ==============================Delete  ===============================================================
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def deleteOrganism(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(Organism, organism_id=pk)
    try:
        object_.delete(**kwargs)
    except Exception as err:
        print(err)
    return redirect('/')
   
# #############################BATCH View############################################
# ---------------------------------------------------------------------------------------------
class BatchCardView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism_Batch 
    template_name = 'dorganism/organism/batch/batch_card.html' 
    filterset_class=Batchfilter
    model_fields=model.HEADER_FIELDS

# ---------------------------------------------------------------------------------------------    
@login_required
def createBatch(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Batch_form()

    if req.method=='POST':
        Organism_Id=req.POST.get('search_organism')
        form=Batch_form(Organism_Id, req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    # print(instance.organism_id)                 
                    instance.save(**kwargs)
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/organism/batch/batch_c.html', { 'form':form, }) 

# ---------------------------------------------------------------------------------------------
from django.http import QueryDict
# ---------------------------------------------------------------------------------------------
@login_required
def updateBatch(req, pk):
    object_=get_object_or_404(Organism_Batch, orgbatch_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Batchupdate_form(instance=object_)
    context={
        "form":form,
        "object":object_,
    }
    if req.method=='PUT':
        qd=QueryDict(req.body).dict()
        print(qd["orgbatch_id"])
        object_batch=object_ 
        form=Batchupdate_form(data=qd, instance=object_batch, )
        
        if form.is_valid():
            kwargs={}
            kwargs['user']=req.user                  
            instance=form.save(commit=False)
            instance.save(**kwargs)
            context={
                "object_batch":object_batch,
                'object':object_batch  # this object refer to the same entry of object_batch
            }
            return render(req, "dorganism/organism/batch/batch_tr.html", context)
        else:
            print(form.errors)
            return render(req, "dorganism/organism/batch/batch_tr.html", context)
    return render(req, "dorganism/organism/batch/batch_u.html", context)

# ---------------------------------------------------------------------------------------------
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def deleteBatch(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(Organism_Batch, orgbatch_id=pk)
    try:
        object_.delete(**kwargs)
    except Exception as err:
        print(err)
    return redirect(req.META['HTTP_REFERER'])  


############################################### Stock View ###########################################

@user_passes_test(lambda u: u.has_permission('Read'), login_url='permission_not_granted') 
def stockList(req, pk):
    # if req.headers.get('x-requested-with') == 'XMLHttpRequest':\
    print(pk)
    res=None
    if req.method == 'GET':
        batch_id=req.GET.get('Batch_id')
        print(f'bactch id {batch_id}')
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
                "stock_date": str(i.stock_date),
                "n_left": str(i.n_left) if i.n_left >1 else " ",
                "n_created": str(i.n_created),
                "stock_notes":str(i.stock_note),
                "biologist": str(i.biologist),
            }
            data.append(item)
        res=data
        
        return JsonResponse({'data':res})
    return JsonResponse({})

# ---------------------------------------------------------------------------------------------
@login_required
def createStock(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Stock_createform()
    if req.method=='POST':
        form=Stock_createform(req.POST)
        if form.is_valid():
            print(req.POST.get("orgbatch_id"))
            orgbatch_id=req.POST.get("orgbatch_id")
            stock_type=req.POST.get("stock_type")
            stock_date=req.POST.get("stock_date")
            kwargs['orgbatch_id']=orgbatch_id
            kwargs['stock_type']=stock_type
            kwargs['stock_date']=stock_date

            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    print("stock saved")
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
        

    return render(req, 'dorganism/organism/batch_stock/stock_c.html', { 'form':form, }) 

# ---------------------------------------------------------------------------------------------
@login_required
def updateStock(req, pk):
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    print(object_)
    kwargs={}
    kwargs['user']=req.user
    form=Stock_form(instance=object_)
    print(pk)
    #-------------------------------------------------------------------------
    if req.method == 'POST' and req.headers.get('x-requested-with') == 'XMLHttpRequest':
        print("start...")
        # process the data sent in the AJAX request
      
        n_left_value=req.POST.get('value')
        object_.n_left=int(n_left_value)-1
        print(object_.n_left)
        
        object_.save(**kwargs)
        print("saved")
        response_data = {'result': 'success'}
        return JsonResponse(response_data)
    if req.method=='POST':
        form=Stock_form(req.POST, instance=object_)
        if "cancel" in req.POST:
            return redirect(req.META['HTTP_REFERER'])
        else:
            try:
                with transaction.atomic(using='dorganism'):        # testing!
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

# ---------------------------------------------------------------------------------------------
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def deleteStock(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    context={'object':object_}
    if req.method=='POST':
        object_.delete(**kwargs)
        return redirect(req.META['HTTP_REFERER'])
    return render(req, "dorganism/organism/batch_stock/stock_d.html", context)

############################################Culture ##################################33
# ==========List View================================Read===========================================
   
@login_required
def createCulture(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Culture_form()

    if req.method=='POST':
        Organism_Id=req.POST.get('search_organism')
        form=Culture_form(Organism_Id, req.POST)
        if form.is_valid():
            print(f'the type is : {req.POST.get("culture_type")}')
            culture_type=req.POST.get("culture_type")
            kwargs['culture_type']=culture_type
            print(f'the source is {req.POST.get("culture_source")}')
            culture_source=req.POST.get("culture_source")
            kwargs['culture_source']=culture_source
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    message=instance.save(**kwargs)
                    
                    if type(message)==Exception:
                        messages.error(req, f'{message} happes')
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/organism/culture/culture_c.html', { 'form':form, }) 

# ---------------------------------------------------------------------------------------------
from django.http import QueryDict

@login_required
def updateCulture(req, pk):
    object_=get_object_or_404(Organism_Culture, id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Cultureupdate_form(instance=object_)
    context={
        "form":form,
        "object":object_,
    }
    if req.method=='PUT':
        qd=QueryDict(req.body).dict()
        object_culture=object_
        form=Cultureupdate_form(data=qd, instance=object_culture )
        
        if form.is_valid():
            kwargs={}
            kwargs['user']=req.user                  
            instance=form.save(commit=False)
            instance.save(**kwargs)
            context={
                "object_cultr":object_culture,
                'object':object_culture  # this object refer to the same entry of object_batch
            }
            return render(req, "dorganism/organism/culture/culture_tr.html", context)
    return render(req, "dorganism/organism/culture/culture_u.html", context)

# ---------------------------------------------------------------------------------------------
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def deleteCulture(req, pk):
    kwargs={}
    kwargs['user']=req.user
    print(f'cultureID {pk}')
    object_=get_object_or_404(Organism_Culture, pk=pk)
    try:
        object_.delete(**kwargs)
        print("delete")
    except Exception as err:
        print(err)
    return redirect(req.META['HTTP_REFERER'])  

