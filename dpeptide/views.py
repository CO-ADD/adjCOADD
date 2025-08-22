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

from applib.django.base.views import (Base_CreateView, Base_UpdateView,  Base_RemoveView, File_CreateView,
                                      Filtered_ListView, permission_not_granted, Htmx_UpdateView)
 
from apputil.utils.upload_steps import UploadHandler_View, SelectSingleFile_StepForm, Upload_StepForm, Finalize_StepForm
from apputil.models import ApplicationLog
from apputil.forms import Document_Form

from adjcoadd.constants import *

from dorganism.models import Taxonomy
from dpeptide.models import  Peptide, Peptide_Batch
from dpeptide.forms import (Peptide_Filter, Peptide_CreateForm, Peptide_UpdateForm)
# , , Cell_UpdateForm, 
#                          CellBatch_Filter, CellBatch_Form, CellBatch_UpdateForm,  
#                          CellBatchStock_Filter, CellBatchStock_Form, CellBatchStock_CreateForm)
#=================================================================================================
# Peptide
#=================================================================================================

class Peptide_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Peptide
    template_name = 'dpeptide/peptide/peptide_list.html'
    filterset_class = Peptide_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Peptide'
    app_name = 'dpeptide'
    ordering=['-acreated_at']
    
# -----------------------------------------------------------------
class Peptide_CardView(Peptide_ListView):
    template_name = 'dpeptide/peptide/peptide_card.html'
    model = Peptide
    model_fields = model.CARDS_FIELDS
    

# -----------------------------------------------------------------
@login_required
#def createCell(req):
def Peptide_CreateView(req):
    '''
    Function View Create new Cell table row with foreignkey: Taxonomy and Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    form=Peptide_CreateForm()
    if req.method=='POST':
        Organism_Name=req.POST.get('search_organism') # -1. Ajax Call(/search_cell/) find Foreignkey cellname
        form=Peptide_CreateForm( Organism_Name, req.POST,) # -2. get create-form with ajax call result
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'): # -3. write new entry in atomic transaction
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Peptide','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dpeptide/peptide/peptide_create.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def Peptide_DetailView(request, pk):
    """
    - Detail view handle Peptide single entry
    display,update and delete.
    - related table overview display, update, create and delete.
    - related table are: batch, stock, culture.
    - data visual table: dataframe and pivot- table
    """
   
    context={}
    object_=get_object_or_404(Peptide, peptide_id=pk)
    try:
        form=Peptide_UpdateForm(initial={'peptide_type':object_.peptide_type, 'peptide_panel':object_.peptide_panel,}, instance=object_)
    except Exception as err:
        print(f"[Peptide_DetailView] {err}")
    context["object"]=object_
    context["form"]=form
    context["doc_form"]=Document_Form
    
    # kwargs={'cell': pk}

    #context["cellbatchimg_form"]=CellBatchImg_Form(cell=pk)

    # data in related tables

    #context["batchimg_obj"]=OrgBatch_Image.objects.filter(orgbatch_id__organism_id=object_.organism_id, astatus__gte=0)
    #context["batchimg_obj_count"]=context["batchimg_obj"].count() if context["batchimg_obj"].count()!=0 else None
   

    # context["batch_obj"]=Cell_Batch.objects.filter(cell_id=object_.cell_id, astatus__gte=0)
    # context["batch_obj_count"]=context["batch_obj"].count() if context["batch_obj"].count()!=0 else None
    # context["batch_fields"]=Cell_Batch.get_fields()

    # context["stock_obj"]=CellBatch_Stock.objects.filter(cellbatch_id__cell_id=object_.cell_id, astatus__gte=0)
    # context["tock_obj_count"]=context["stock_obj"].count() if context["stock_obj"].count()!=0 else None
    # context["tock_fields"]=CellBatch_Stock.get_fields()

    # context["cell_stock_count"]=Cell_Batch.objects.annotate(number_of_stocks=Count('cellbatch_id')) 
    #context["cultr_obj"]=Cell_Culture.objects.filter(cell_id=object_.cell_id, astatus__gte=0)
    #context["cultr_obj_count"]=context["cultr_obj"].count() if context["cultr_obj"].count()!=0 else None
    #context["cultr_fields"]=Cell_Culture.get_fields() 
    # if 'cell_id' in context["cultr_fields"]:
    #     context["cultr_fields"].remove('cell_id')    # customize LIST_VIEW_FIELDS
    # context["vitekast_obj"]=SimpleLazyObject(lambda: VITEK_AST.objects.filter(organism=object_.organism_name, astatus__gte=0))
    # context["vitekast_obj_count"]=context["vitekast_obj"].count() if context["vitekast_obj"].count()!=0 else None
    # context["vitekast_fields"]=VITEK_AST.get_fields(fields=VITEK_AST.LIST_VIEW_FIELDS)

    # data in pivotted and highlighted Tables
    
    # if request.method == 'POST':
    #     displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source'] #<what is needed from here?>
    #     context["table"] = data_frame_style(pk, displaycols)['style_table']
    #     context["df_entries"] = data_frame_style(pk, displaycols)['df_entries']
    #     context["pivottable"] = pivottable_style(pk)
    #     return render(request, "dcell/cell/cell_mic.html", context)
    
    return render(request, "dpeptide/peptide/peptide_detail.html", context)

# -----------------------------------------------------------------
@login_required
def Peptide_UpdateView(req, pk):
    object_=get_object_or_404(Peptide, peptide_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Peptide_UpdateForm(initial={'peptide_type':object_.peptide_type, 'peptide_panel':object_.peptide_panel, 'assoc_documents': [i.doc_file for i in object_.assoc_documents.all()]}, instance=object_)
    if object_.organism_name.org_class: # Organism_Class_str for display class
        Organism_Class_str=object_.organism_name.org_class.dict_value
    else:
        Organism_Class_str="No Class"
    #print(f"[Cell_UpdateView] Organism_Class_str")
    if req.method=='POST':
        #print(f"[Cell_UpdateView] POST")
        try:
            with transaction.atomic(using='dpeptide'):        # testing!
                obj = Peptide.objects.select_for_update().get(peptide_id=pk)
                #If update Cell Name
                if  req.POST.get('search_peptide'):
                    Organism_Name_str=req.POST.get('search_peptide')
                    Organism_new_obj=get_object_or_404(Taxonomy, organism_name=Organism_Name_str)
                    
                    form=Peptide_UpdateForm(Organism_Name_str, req.POST, instance=obj)
                    #-Not allow to update name in different class--
                    if Organism_new_obj.org_class.dict_value and Organism_new_obj.org_class.dict_value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                else:
                    #print("Else")
                    form=Peptide_UpdateForm(object_.organism_name, req.POST, instance=obj) 
                
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
   
    return render(req, "dpeptide/peptide/peptide_update.html", context)

# -----------------------------------------------------------------
class Peptide_RemoveView(Base_RemoveView):
    model = Peptide
    transaction_use = 'ddpeptide'