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
from django.db.models import Count
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.utils.functional import SimpleLazyObject

from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from apputil.utils.filters_base import FilteredListView
from apputil.utils.views_base import permission_not_granted, HtmxupdateView, SimplecreateView, SimpleupdateView,  SimpledeleteView, CreateFileView

from adjcoadd.constants import *

from dorganism.models import  Taxonomy, Organism, Organism_Batch, OrgBatch_Stock, Organism_Culture, OrgBatch_Image
from dorganism.forms import (Taxonomy_Filter, Taxonomy_Form,
                            Organism_Filter, CreateOrganism_form, UpdateOrganism_form, 
                            OrgBatch_Filter, OrgBatch_Form, OrgBatch_UpdateForm,  
                            OrgBatchStock_Form, OrgBatchStock_Filter, OrgBatchStock_CreateForm, 
                            OrgCulture_Form, OrgCulture_UpdateForm,
                            OrgBatchImg_Form,)
from ddrug.models import VITEK_AST, MIC_COADD
from ddrug.utils.antibiogram import get_Antibiogram_byOrgID_Html
#from dorganism.utils.data_visual import data_frame_style, pivottable_style
    
#=================================================================================================
# Taxonomy
#=================================================================================================

class Taxonomy_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Taxonomy  
    template_name = 'dorganism/taxonomy/taxonomy_list.html' 
    filterset_class=Taxonomy_Filter
    model_fields=model.HEADER_FIELDS
    model_name = 'Taxonomy'
    app_name = 'dorganism'
    ordering = ['org_class']

# -----------------------------------------------------------------
class Taxonomy_CardView(Taxonomy_ListView):
    template_name = 'dorganism/taxonomy/taxonomy_card.html'
    
# -----------------------------------------------------------------
@login_required
def Taxonomy_DetailView(req, slug=None):
    context={}
    object_=get_object_or_404(Taxonomy, urlname=slug)
    context["object"]=object_
    context['form']=Taxonomy_Form(instance=object_)
    return render(req, "dorganism/taxonomy/taxonomy_detail.html", context)

# -----------------------------------------------------------------
class Taxonomy_CreateView(SimplecreateView):
    form_class=Taxonomy_Form
    template_name='dorganism/taxonomy/taxonomy_c.html'
    transaction_use = 'dorganism'
        
# -----------------------------------------------------------------
class Taxonomy_UpdateView(SimpleupdateView):
    form_class=Taxonomy_Form
    template_name='dorganism/taxonomy/taxonomy_u.html'
    model=Taxonomy
    transaction_use = 'dorganism'

# -----------------------------------------------------------------
class Taxonomy_DeleteView(SimpledeleteView):
    model = Taxonomy
    transaction_use = 'dorganism'

#=================================================================================================
# Organism
#=================================================================================================

class Organism_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = Organism  
    template_name = 'dorganism/organism/organism_list.html'
    filterset_class = Organism_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'Organism'
    app_name = 'dorganism'
    ordering=['-acreated_at']
    
# -----------------------------------------------------------------
class Organism_CardView(Organism_ListView):
    template_name = 'dorganism/organism/organism_card.html'
    model = Organism  
    model_fields = model.CARDS_FIELDS


# -----------------------------------------------------------------
@login_required
def Organism_CreateView(req):
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
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Organism','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dorganism/organism/organism_c.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def Organism_DetailView(request, pk):
    """
    - Detail view handle organism single entry
    display,update and delete.
    - related table overview display, update, create and delete.
    - related table are: batch, stock, culture.
    - data visual table: dataframe and pivot- table
    """
   
    context={}
    object_=get_object_or_404(Organism, organism_id=pk)
    try:
        form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel,}, instance=object_)
    except Exception as err:
        print(err)
    context["object"]=object_
    context["form"]=form
    context["doc_form"]=Document_Form
    # kwargs={'org': pk}

    context["orgbatchimg_form"]=OrgBatchImg_Form(org=pk)

    # data in related tables

    context["batchimg_obj"]=OrgBatch_Image.objects.filter(orgbatch_id__organism_id=object_.organism_id, astatus__gte=0)
    context["batchimg_obj_count"]=context["batchimg_obj"].count() if context["batchimg_obj"].count()!=0 else None
   

    context["batch_obj"]=Organism_Batch.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["batch_obj_count"]=context["batch_obj"].count() if context["batch_obj"].count()!=0 else None
    context["batch_fields"]=Organism_Batch.get_fields()

    # context["stock_obj"]=OrgBatch_Stock.objects.filter(orgbatch_id__organism_id=object_.organism_id, astatus__gte=0)
    # context["stock_obj_count"]=context["stock_obj"].count() if context["stock_obj"].count()!=0 else None
    # context["stock_fields"]=OrgBatch_Stock.get_fields()

    context["stock_count"]=Organism_Batch.objects.annotate(number_of_stocks=Count('orgbatch_id')) 
    context["cultr_obj"]=Organism_Culture.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["cultr_obj_count"]=context["cultr_obj"].count() if context["cultr_obj"].count()!=0 else None
    context["cultr_fields"]=Organism_Culture.get_fields() 
    if 'organism_id' in context["cultr_fields"]:
        context["cultr_fields"].remove('organism_id')    # customize HEADER_FIELDS
    context["vitekast_obj"]=SimpleLazyObject(lambda: VITEK_AST.objects.filter(organism=object_.organism_name, astatus__gte=0))
    context["vitekast_obj_count"]=context["vitekast_obj"].count() if context["vitekast_obj"].count()!=0 else None
    context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEK_AST.HEADER_FIELDS)
    context["n_entries"] = 0

    # data in pivotted and highlighted Tables
    if request.method == 'POST':
        displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source']
        print('POST')
        _pivDict = get_Antibiogram_byOrgID_Html(pk, displaycols)
        print(_pivDict['n_entries'])
        context["table"] = _pivDict['html_table']
        context["n_entries"] = _pivDict['n_entries']
        context["pivottable"] = _pivDict['pivot_table']
        return render(request, "dorganism/organism/organism_mic.html", context)
    
    return render(request, "dorganism/organism/organism_detail.html", context)

# -----------------------------------------------------------------
@login_required
def Organism_UpdateView(req, pk):
    object_=get_object_or_404(Organism, organism_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=UpdateOrganism_form(initial={'strain_type':object_.strain_type, 'strain_panel':object_.strain_panel, 'assoc_documents': [i.doc_file for i in object_.assoc_documents.all()]}, instance=object_)
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
                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update Organism','Completed')
                    # form.save_m2m() 
                    return redirect(req.META['HTTP_REFERER'])
                else:
                    messages.warning(req, f'Update failed due to {form.errors} error')
                   
        except Exception as err:
            print(err)
            messages.warning(req, f'Update failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])
  
    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dorganism/organism/organism_u.html", context)

# -----------------------------------------------------------------
class Organism_DeleteView(SimpledeleteView):
    model = Organism
    transaction_use = 'dorganism'

#=================================================================================================
# OrgBatch  
#=================================================================================================

class OrgBatch_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism_Batch 
    template_name = 'dorganism/orgbatch/orgbatch_list.html' 
    filterset_class=OrgBatch_Filter
    model_fields=model.HEADER_FIELDS
    model_name = 'Organism_Batch'
    app_name = 'dorganism'

# -----------------------------------------------------------------
@login_required
def OrgBatch_CreateView(req, organism_id):
    kwargs={'user': req.user}
    form=OrgBatch_Form()
    
    if req.method=='POST':
        form=OrgBatch_Form(req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.organism_id=get_object_or_404(Organism, pk=organism_id)              
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new OrgBatch','Completed')
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/orgbatch/orgbatch_c.html', { 'form':form, 'organism_id':organism_id}) 

# -----------------------------------------------------------------
class OrgBatch_UpdateView(HtmxupdateView):
    form_class=OrgBatch_UpdateForm
    template_name="dorganism/orgbatch/orgbatch_u.html"
    template_partial="dorganism/orgbatch/orgbatch_tr.html"
    model=Organism_Batch
    transaction_use = 'dorganism'

# -----------------------------------------------------------------
class OrgBatch_DeleteView(SimpledeleteView):
    model = Organism_Batch
    transaction_use = 'dorganism'

#=================================================================================================
# OrgBatch Stock  
#=================================================================================================
class OrgBatchStock_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = OrgBatch_Stock  
    template_name = 'dorganism/orgbatchstock/orgbatchstock_list.html'
    filterset_class = OrgBatchStock_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'OrgBatch_Stock'
    app_name = 'dorganism'



#-------------------------------------------------------------------------------
# within Organism_DetailView
## here is response to an Ajax call
## to send data to child datatable 
@user_passes_test(lambda u: u.has_permission('Read'), login_url='permission_not_granted') 
def OrgBatchStock_DetailView(req, pk):
    res=None
    if req.method == 'GET':
        batch_id=req.GET.get('Batch_id')
        object_=get_object_or_404(Organism_Batch, orgbatch_id=batch_id)#Organism_Batch.objects.get(orgbatch_id=batch_id)
        print(object_)
        qs=OrgBatch_Stock.objects.filter(orgbatch_id=object_, astatus__gte=0, n_left__gt=0) # n_left show when bigger or equal to 2
        data=[]
        for i in qs:
            item={
                "stock_id":i.pk,
                "stock_type": str(i.stock_type.dict_value) or 'no data',
                "location_freezer":str(i.location_freezer) or 'no data',
                "location_rack": str(i.location_rack),
                "location_col": str(i.location_column),
                "location_slot": str(i.location_slot),
                "stock_date": str(i.stock_date.strftime("%d-%m-%Y")) if i.stock_date else '-',
                "n_left": str(i.n_left),
                "n_created": str(i.n_created),
                "stock_notes":str(i.stock_note),
                "biologist": str(i.biologist),
            }
            data.append(item)          
        res=data        
        return JsonResponse({'data':res})
    return JsonResponse({})
      
#-------------------------------------------------------------------------------
@login_required
def OrgBatchStock_CreateView(req, orgbatch_id):
    kwargs={}
    kwargs['user']=req.user
    form = OrgBatchStock_CreateForm(initial={"orgbatch_id":orgbatch_id},)
    if req.method=='POST':
        form=OrgBatchStock_CreateForm(req.POST)
        if form.is_valid():

            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create new OrgBatchStock','Completed')
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'wrong {form.errors}')
    return render(req, 'dorganism/orgbatchstock/orgbatchstock_c.html', { 'form':form, 'orgbatch_id':orgbatch_id }) 

#-------------------------------------------------------------------------------
@login_required
def OrgBatchStock_UpdateView(req, pk):
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    kwargs={}
    kwargs['user']=req.user
    form=OrgBatchStock_Form(instance=object_)
    if req.method == 'POST' and req.headers.get('x-requested-with') == 'XMLHttpRequest':
        # process the data sent in the AJAX request
        n_left_value=req.POST.get('value')
        object_.n_left=int(n_left_value)-1
        object_.save(**kwargs)
        ApplicationLog.add('Updated',str(object_.pk),'Info',req.user,str(object_.pk),'Updated Stock_N_Left','Completed')
        response_data = {'result': str(object_.n_left)}
        return JsonResponse(response_data)
    if req.method=='POST':
        form=OrgBatchStock_Form(req.POST, instance=object_)
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
                            ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update OrgBatchStock','Completed')
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
    return render(req, "dorganism/orgbatchstock/orgbatchstock_u.html", context)

#-------------------------------------------------------------------------------
class OrgBatchStock_DeleteView(SimpledeleteView):
    model = OrgBatch_Stock
    transaction_use = 'dorganism'

#=================================================================================================
# Organim Culture 
#=================================================================================================

@login_required
def OrgCulture_CreateView(req, organism_id):
    
    kwargs={'user' : req.user}
    form=OrgCulture_Form()
    
    if req.method=='POST':
        form=OrgCulture_Form(req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False)
                    instance.organism_id=get_object_or_404(Organism, pk=organism_id)
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create new OrgCulture','Completed')                  
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/orgculture/orgculture_c.html', { 'form':form, 'organism_id':organism_id}) 

#-------------------------------------------------------------------------------
class OrgCulture_UpdateView(HtmxupdateView):
    form_class=OrgCulture_UpdateForm
    template_name="dorganism/orgculture/orgculture_u.html"
    template_partial="dorganism/orgculture/orgculture_tr.html"
    model=Organism_Culture
    transaction_use = 'dorganism'

#-------------------------------------------------------------------------------
class OrgCulture_DeleteView(SimpledeleteView):
    model = Organism_Culture
    transaction_use = 'dorganism'

#=================================================================================================
# OrgBatchImage OrgBatchImage OrgBatchImage OrgBatchImage OrgBatchImage OrgBatchImage OrgBatchImage
#=================================================================================================
class OrgBatchImg_DeleteView(SimpledeleteView):
    model = OrgBatch_Image
    transaction_use = 'dorganism'

#-------------------------------------------------------------------------------
class OrgBatchImg_CreateView(CreateFileView):
    form_class=OrgBatchImg_Form
    model = Organism
    file_field = 'image_file' #this is uploading field name of Orgbatchimg
    transaction_use = 'dorganism'

    def form_valid(self, form):
        instance = form.save(commit=False)
        if getattr(instance, self.file_field):
            with transaction.atomic(using=self.transaction_use):
                kwargs={'user': self.request.user}
                instance.save(**kwargs)
        else:
            messages.warning(self.request, 'No file provided.')
        return redirect(self.request.META['HTTP_REFERER'])
    
    # def get_form_kwargs(self):
    #     kwargs = super().get_form_kwargs()
    #     print("1",self.request)
    #     kwargs["pk"] = self.request
    #     return kwargs

    # def get_form(self, form_class=None):
    #     print("wil called")
    #     form = super().get_form(form_class)
    #     organism = get_object_or_404(self.model, organism_id=kwargs['pk'])
    #     form.fields['orgbatch_id'].queryset = Organism_Batch.objects.filter(organism_id = organism.pk)
    #     return form     
