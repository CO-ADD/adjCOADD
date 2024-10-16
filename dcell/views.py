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
from apputil.utils.upload_steps import UploadHandler_View, SelectSingleFile_StepForm, Upload_StepForm, Finalize_StepForm

from adjcoadd.constants import *

from dorganism.models import Taxonomy
from dcell.models import  Cell, Cell_Batch, CellBatch_Stock
from dcell.forms import (Cell_Filter, Cell_CreateForm, Cell_UpdateForm, 
                         CellBatch_Filter, CellBatch_Form, CellBatch_UpdateForm,  
                         CellBatchStock_Filter, CellBatchStock_Form, CellBatchStock_CreateForm)

from dcell.utils.upload_cell import upload_Cells_Process

#from ddrug.models import VITEK_AST, MIC_COADD
#from dorganism.utils.data_visual import data_frame_style, pivottable_style
    

#=================================================================================================
# Cell
#=================================================================================================

class Cell_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = Cell
    template_name = 'dcell/cell/cell_list.html'
    filterset_class = Cell_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'Cell'
    app_name = 'dcell'
    ordering=['-acreated_at']
    
# -----------------------------------------------------------------
class Cell_CardView(Cell_ListView):
    template_name = 'dcell/cell/cell_card.html'
    model = Cell
    model_fields = model.CARDS_FIELDS


# -----------------------------------------------------------------
@login_required
#def createCell(req):
def Cell_CreateView(req):
    '''
    Function View Create new Cell table row with foreignkey: Taxonomy and Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    form=Cell_CreateForm()
    if req.method=='POST':
        Organism_Name=req.POST.get('search_organism') # -1. Ajax Call(/search_cell/) find Foreignkey cellname
        form=Cell_CreateForm( Organism_Name, req.POST,) # -2. get create-form with ajax call result
        if form.is_valid():
            try:
                with transaction.atomic(using='dcell'): # -3. write new entry in atomic transaction
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Cell','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dcell/cell/cell_c.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def Cell_DetailView(request, pk):
    """
    - Detail view handle cell single entry
    display,update and delete.
    - related table overview display, update, create and delete.
    - related table are: batch, stock, culture.
    - data visual table: dataframe and pivot- table
    """
   
    context={}
    object_=get_object_or_404(Cell, cell_id=pk)
    try:
        form=Cell_UpdateForm(initial={'cell_type':object_.cell_type, 'cell_panel':object_.cell_panel,}, instance=object_)
    except Exception as err:
        print(f"[Cell_DetailView] {err}")
    context["object"]=object_
    context["form"]=form
    context["doc_form"]=Document_Form
    # kwargs={'cell': pk}

    #context["cellbatchimg_form"]=CellBatchImg_Form(cell=pk)

    # data in related tables

    #context["batchimg_obj"]=OrgBatch_Image.objects.filter(orgbatch_id__organism_id=object_.organism_id, astatus__gte=0)
    #context["batchimg_obj_count"]=context["batchimg_obj"].count() if context["batchimg_obj"].count()!=0 else None
   

    context["batch_obj"]=Cell_Batch.objects.filter(cell_id=object_.cell_id, astatus__gte=0)
    context["batch_obj_count"]=context["batch_obj"].count() if context["batch_obj"].count()!=0 else None
    context["batch_fields"]=Cell_Batch.get_fields()

    context["stock_obj"]=CellBatch_Stock.objects.filter(cellbatch_id__cell_id=object_.cell_id, astatus__gte=0)
    context["tock_obj_count"]=context["stock_obj"].count() if context["stock_obj"].count()!=0 else None
    context["tock_fields"]=CellBatch_Stock.get_fields()

    context["cell_stock_count"]=Cell_Batch.objects.annotate(number_of_stocks=Count('cellbatch_id')) 
    #context["cultr_obj"]=Cell_Culture.objects.filter(cell_id=object_.cell_id, astatus__gte=0)
    #context["cultr_obj_count"]=context["cultr_obj"].count() if context["cultr_obj"].count()!=0 else None
    #context["cultr_fields"]=Cell_Culture.get_fields() 
    # if 'cell_id' in context["cultr_fields"]:
    #     context["cultr_fields"].remove('cell_id')    # customize HEADER_FIELDS
    # context["vitekast_obj"]=SimpleLazyObject(lambda: VITEK_AST.objects.filter(organism=object_.organism_name, astatus__gte=0))
    # context["vitekast_obj_count"]=context["vitekast_obj"].count() if context["vitekast_obj"].count()!=0 else None
    # context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEK_AST.HEADER_FIELDS)

    # data in pivotted and highlighted Tables
    
    # if request.method == 'POST':
    #     displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source'] #<what is needed from here?>
    #     context["table"] = data_frame_style(pk, displaycols)['style_table']
    #     context["df_entries"] = data_frame_style(pk, displaycols)['df_entries']
    #     context["pivottable"] = pivottable_style(pk)
    #     return render(request, "dcell/cell/cell_mic.html", context)
    
    return render(request, "dcell/cell/cell_detail.html", context)

# -----------------------------------------------------------------
@login_required
def Cell_UpdateView(req, pk):
    object_=get_object_or_404(Cell, cell_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Cell_UpdateForm(initial={'cell_type':object_.cell_type, 'cell_panel':object_.cell_panel, 'assoc_documents': [i.doc_file for i in object_.assoc_documents.all()]}, instance=object_)
    if object_.organism_name.org_class: # Organism_Class_str for display class
        Organism_Class_str=object_.organism_name.org_class.dict_value
    else:
        Organism_Class_str="No Class"
    #print(f"[Cell_UpdateView] Organism_Class_str")
    if req.method=='POST':
        #print(f"[Cell_UpdateView] POST")
        try:
            with transaction.atomic(using='dcell'):        # testing!
                obj = Cell.objects.select_for_update().get(cell_id=pk)
                #If update Cell Name
                if  req.POST.get('search_cell'):
                    Organism_Name_str=req.POST.get('search_cell')
                    Organism_new_obj=get_object_or_404(Taxonomy, organism_name=Organism_Name_str)
                    
                    form=Cell_UpdateForm(Organism_Name_str, req.POST, instance=obj)
                    #-Not allow to update name in different class--
                    if Organism_new_obj.org_class.dict_value and Organism_new_obj.org_class.dict_value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                else:
                    #print("Else")
                    form=Cell_UpdateForm(object_.organism_name, req.POST, instance=obj) 
                
                if form.is_valid():  
                    #print(f"Saving {obj}")     
                    instance=form.save(commit=False)
                    instance.save(**kwargs)
                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Updated Cell','Completed')
                    # form.save_m2m() 
                    return redirect(req.META['HTTP_REFERER'])
                else:
                
                    messages.warning(req, f'Update 1 failed due to {form.errors} error')
                   
        except Exception as err:
            messages.warning(req, f'Update 2 failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])
  
    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dcell/cell/cell_u.html", context)

# -----------------------------------------------------------------
class Cell_DeleteView(SimpledeleteView):
    model = Cell
    transaction_use = 'dcell'


# -----------------------------------------------------------------
class Cell_Upload_HandlerView(UploadHandler_View):
    
    name_step1="Upload"
    form_list = [
        ('select_file', SelectSingleFile_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dcell/cell/importhandler_cell.html'
    
    # customize util functions to validate files:
    # vitek -- upload_VitekPDF_Process
    def file_process_handler(self, request, *args, **kwargs):
        try:
            form_data=kwargs.get('form_data', None)
        except Exception as err:
            print(f"[Import_CellView] {err}")
        
        valLog=upload_Cells_Process(request, self.dirname, self.filelist, upload=self.upload, appuser=request.user) 
        return(valLog)

#=================================================================================================
# CellBatch  
#=================================================================================================

class CellBatch_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Cell_Batch 
    template_name = 'dcell/cellbatch/cellbatch_list.html' 
    filterset_class=CellBatch_Filter
    model_fields=model.HEADER_FIELDS
    model_name = 'Cell_Batch'
    app_name = 'dcell'

# -----------------------------------------------------------------
@login_required
def CellBatch_CreateView(req, cell_id):
    kwargs={'user': req.user}
    form=CellBatch_Form()
    
    if req.method=='POST':
        form=CellBatch_Form(req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dcell'):
                    instance=form.save(commit=False) 
                    instance.cell_id=get_object_or_404(Cell, pk=cell_id)              
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a CellBatch','Completed')
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dcell/cellbatch/cellbatch_c.html', { 'form':form, 'cell_id':cell_id}) 

# -----------------------------------------------------------------
class CellBatch_UpdateView(HtmxupdateView):
    form_class=CellBatch_UpdateForm
    template_name="dcell/cellbatch/cellbatch_u.html"
    template_partial="dcell/cellbatch/cellbatch_tr.html"
    model=Cell_Batch
    transaction_use = 'dcell'

# -----------------------------------------------------------------
class CellBatch_DeleteView(SimpledeleteView):
    model = Cell_Batch
    transaction_use = 'dcell'

#=================================================================================================
# CellBatch Stock  
#=================================================================================================
class CellBatchStock_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = CellBatch_Stock  
    template_name = 'dcell/cellbatchstock/cellbatchstock_list.html'
    filterset_class = CellBatchStock_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'CellBatch_Stock'
    app_name = 'dcell'



#-------------------------------------------------------------------------------
# within Cell_DetailView
## here is response to an Ajax call
## to send data to child datatable 
@user_passes_test(lambda u: u.has_permission('Read'), login_url='permission_not_granted') 
def CellBatchStock_DetailView(req, pk):
    res=None
    if req.method == 'GET':
        batch_id=req.GET.get('Batch_id')
        object_=get_object_or_404(Cell_Batch, cellbatch_id=batch_id)#Cell_Batch.objects.get(cellbatch_id=batch_id)
        #print(f"[CellBatchStock_DetailView] {object_}"")
        qs=CellBatch_Stock.objects.filter(cellbatch_id=object_, astatus__gte=0, n_left__gt=0) # n_left show when bigger or equal to 2
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
def CellBatchStock_CreateView(req, cellbatch_id):
    kwargs={}
    kwargs['user']=req.user
    form = CellBatchStock_CreateForm(initial={"cellbatch_id":cellbatch_id},)
    if req.method=='POST':
        form=CellBatchStock_CreateForm(req.POST)
        if form.is_valid():

            try:
                with transaction.atomic(using='dcell'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new CellBatchStock','Completed')
                    return redirect(req.META['HTTP_REFERER']) 
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'[CellBatchStock_CreateView] {form.errors}')
    return render(req, 'dcell/cellbatchstock/cellbatchstock_c.html', { 'form':form, 'cellbatch_id':cellbatch_id }) 

#-------------------------------------------------------------------------------
@login_required
def CellBatchStock_UpdateView(req, pk):
    object_=get_object_or_404(CellBatch_Stock, pk=pk)
    kwargs={}
    kwargs['user']=req.user
    form=CellBatchStock_Form(instance=object_)
    if req.method == 'POST' and req.headers.get('x-requested-with') == 'XMLHttpRequest':
        # process the data sent in the AJAX request
        n_left_value=req.POST.get('value')
        object_.n_left=int(n_left_value)-1
        object_.save(**kwargs)
        ApplicationLog.add('Updated',str(object_.pk),'Info',req.user,str(object_.pk),'Updated Stock_N_Left','Completed')
        response_data = {'result': str(object_.n_left)}
        return JsonResponse(response_data)
    if req.method=='POST':
        form=CellBatchStock_Form(req.POST, instance=object_)
        if "cancel" in req.POST:
            return redirect(req.META['HTTP_REFERER'])
        else:
            try:
                with transaction.atomic(using='dcell'):      
                    obj = CellBatch_Stock.objects.select_for_update().get(pk=pk)
                    try:
                        if form.is_valid():               
                            instance=form.save(commit=False)
                            instance.save(**kwargs)
                            ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Updated CellBatchStock','Completed')
                            return redirect(req.META['HTTP_REFERER'])
                            
                    except Exception as err:
                        print(f'[CellBatchStock_UpdateView] form.error:  {form.errors}, error: {err}')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
                return redirect(req.META['HTTP_REFERER'])
    context={
        "form":form,
        "object":object_,
    }
    return render(req, "dcell/cellbatchstock/cellbatchstock_u.html", context)

#-------------------------------------------------------------------------------
class CellBatchStock_DeleteView(SimpledeleteView):
    model = CellBatch_Stock
    transaction_use = 'dcell'
