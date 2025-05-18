
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
from applib.django.views import Base_CreateView, Base_UpdateView, Base_DeleteView, Filtered_ListView
from apputil.utils.form_wizard_tools import ImportHandler_View, SelectMultipleFiles_StepForm, Upload_StepForm, Finalize_StepForm

# from apputil.utils.filters_base import FilteredListView
# from apputil.utils.views_base import permission_not_granted, HtmxupdateView, SimplecreateView, SimpleupdateView,  SimpledeleteView, CreateFileView

# from adjcoadd.constants import *

from dscreen.models import Screen_Run
from dscreen.forms import ScreenRun_Filter, ScreenRun_CreateForm, ScreenRun_UpdateForm
from dsample.models import Project
from dplate.models import MasterPlate, TestPlate

#=================================================================================================
# ScreenRun
#=================================================================================================
class ScreenRun_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Screen_Run  
    template_name = 'dscreen/screenrun/screenrun_list.html'
    filterset_class = ScreenRun_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'Screen_Run'
    app_name = 'dscreen'
    ordering=['-acreated_at']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = 'coadd_base.html'
        return context

# -----------------------------------------------------------------
# class ScreenRun_CardView(ScreenRun_ListView):
#     template_name = 'dscreen/screenrun/screenrun_card.html'
#     model = Screen_Run  
#     model_fields = model.CARDS_FIELDS

# -----------------------------------------------------------------
@login_required
def ScreenRun_CreateView(req):
    '''
    View to Create new ScreenRun foreignkey: Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    form=ScreenRun_CreateForm()
    if req.method=='POST':
        form=ScreenRun_CreateForm(req.POST) 
        if form.is_valid():
            try:
                with transaction.atomic(using='dscreen'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Screen Run','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dscreen/screenrun/screenrun_create.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def ScreenRun_DetailView(req, pk):
    """
    - Detail view handle ScreenRun entry display,update and delete.
    - related table overview display.
    - related table are: testplate, masterplate, Processing.
    - data visual table: dataframe and pivot- table
    """
    context={}
    # try:
    _object=get_object_or_404(Screen_Run, run_id=pk)
    form=ScreenRun_UpdateForm(initial={'run_type':_object.run_type, 
                                      'run_status':_object.run_status,}, 
                                    instance=_object)
    print(f"[ScreenRun_DetailView] {req.method}")
    if req.method == 'POST':
        print(f"[ScreenRun_DetailView] {req.POST}")

    context["object"]=_object
    context["form"]=form

    # plate_data_df = get_screenrun_plates(_object.run_id)
    # context["org_id_obj_count"] = len(id_data_df)
    # context["org_id_obj"] = id_data_df.values.tolist()
    # context["org_id_fields"] = list(id_data_df.columns)

    # project_data_df = get_screenrun_projects(_object.run_id)
    # context["org_id_obj_count"] = len(id_data_df)
    # context["org_id_obj"] = id_data_df.values.tolist()
    # context["org_id_fields"] = list(id_data_df.columns)

    return render(req, "dscreen/screenrun/screenrun_detail.html", context)

# -----------------------------------------------------------------
@login_required
def ScreenRun_UpdateView(req, pk):
    _object=get_object_or_404(Screen_Run, run_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=ScreenRun_UpdateForm(initial={'run_type':_object.run_type, 
                                      'run_status':_object.run_status,}, 
                                    instance=_object)
    print(f"[ScreenRun_UpdateView] {req.method}")
    print(f"[ScreenRun_UpdateView] {req.session}")
    if req.method=='POST':
        print(f"[ScreenRun_UpdateView] {req.POST}")
        try:
            with transaction.atomic(using='dscreen'):
                obj = Screen_Run.objects.select_for_update().get(run_id=pk)
                form=ScreenRun_UpdateForm(req.POST, instance=obj)    
                if form.is_valid():       
                    instance=form.save(commit=False)
                    instance.save(**kwargs)
                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update Screen_Run','Completed')
                    # form.save_m2m() 
                    return redirect(req.META['HTTP_REFERER'])
                else:
                    messages.warning(req, f'Update failed due to {form.errors} error')
                    
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])

    context={}
    context["object"]=_object
    context["form"]=form
   
    return render(req, "dscreen/screenrun/screenrun_update.html", context)

# -----------------------------------------------------------------
class ScreenRun_DeleteView(Base_DeleteView):
    model = Screen_Run
    transaction_use = 'dscreen'

@login_required
def Add_Readouts(req, pk):
    context = {}
    _object = get_object_or_404(Screen_Run, run_id=pk)

    form=ScreenRun_UpdateForm(instance=_object)
    context["object"]=_object
    context["form"] = form
    return render(req,'dscreen/screenrun/add_readouts.html',context)


# -----------------------------------------------------------------
# class Add_Readouts(ImportHandler_View):
#     model = Screen_Run

#     name_step1="Upload"
#     form_list = [
#         ('select_file', SelectMultipleFiles_StepForm),
#         #('upload', VitekValidation_StepForm),
#         ('finalize', Finalize_StepForm),
#     ]
#     template_name = 'dscreen/screenrun_add_readouts.html'


    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.run_id=None
    
    # # customize util functions to validate files:
    # # vitek -- upload_VitekPDF_Process
    # def file_process_handler(self, request, *args, **kwargs):
    #     try:
    #         form_data=kwargs.get('form_data', None)
    #     except Exception as err:
    #         print(err)
    #         return (err)
    #     if 'upload-orgbatch_id' in form_data.keys():
    #         self.organism_batch=form_data['upload-orgbatch_id'] #get organism_batch  
    #         print(self.organism_batch)   
    #     valLog=upload_VitekPDF_Process(request, self.dirname, self.filelist, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user) 
    #     return(valLog)
