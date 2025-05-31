import os
import json
import datetime
from io import BytesIO as IO

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
from django.urls import reverse_lazy, reverse
from django.utils.functional import SimpleLazyObject

from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_DeleteView, Filtered_ListView
 
# from applib.django.base.filters import Filtered_ListView
# from applib.django.base.views import permission_not_granted, Htmx_UpdateView, Base_CreateView, Base_UpdateView,  Base_DeleteView, File_CreateView

# from adjcoadd.constants import *

from dsample.models import Project
from dsample.forms import Project_Filter, Project_CreateForm, Project_UpdateForm
from dscreen.models import Screen_Run
from dplate.models import MasterPlate, TestPlate
from applib.report.screen_data import Report_Screening


#=================================================================================================
# Project
#=================================================================================================
class Project_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Project  
    template_name = 'dsample/project/project_list.html'
    filterset_class = Project_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'Project'
    app_name = 'dsample'
    ordering=['-acreated_at']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = 'coadd_base.html'
        return context

# -----------------------------------------------------------------
@login_required
def Project_CreateView(req):
    '''
    View to Create new Project foreignkey: Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    form=Project_CreateForm()
    if req.method=='POST':
        form=Project_CreateForm(req.POST) 
        if form.is_valid():
            print('Project_CreateView Valid')
            try:
                with transaction.atomic(using='dsample'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Project','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dsample/project/project_create.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def Project_DetailView(req, pk):
    """
    - Detail view handle Project entry display,update and delete.
    - related table overview display.
    - related table are: testplate, masterplate, Processing.
    - data visual table: dataframe and pivot- table
    """
    context={}
    # try:
    _object=get_object_or_404(Project, project_id=pk)
    form=Project_UpdateForm(instance=_object)
    # form=Project_UpdateForm(initial={'project_type':_object.project_type, 
    #                                   'project_status':_object.project_status,}, 
    #                                 instance=_object)
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

    return render(req, "dsample/project/project_detail.html", context)

# -----------------------------------------------------------------
@login_required
def Project_UpdateView(req, pk):
    _object=get_object_or_404(Project, project_id=pk)
    kwargs={}
    kwargs['user']=req.user
    form=Project_UpdateForm(initial={'project_type':_object.project_type, 
                                      'project_status':_object.project_status,}, 
                                    instance=_object)
    if req.method=='POST':
        try:
            with transaction.atomic(using='dsample'):
                obj = Project.objects.select_for_update().get(project_id=pk)
                form=Project_UpdateForm(req.POST, instance=obj)    
                if form.is_valid():       
                    instance=form.save(commit=False)
                    instance.save(**kwargs)
                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update Project','Completed')
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
   
    return render(req, "dsample/project/project_update.html", context)

# -----------------------------------------------------------------
class Project_DeleteView(Base_DeleteView):
    model = Project
    transaction_use = 'dsample'


# -----------------------------------------------------------------
@login_required
def Project_ReportView(req, pk):
    

    if req.method=='GET':

        _object=get_object_or_404(Project, project_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'Project_{pk}_Summary_{_now:%Y%m%d}.xlsx'

        cReport = Report_Screening()
        cReport.qry_by_ProjectID(_object)
        if cReport.n_compounds>0:
            cReport.get_dataframe()
            cReport.get_sample_info(Storage_Info=False, Structure_Info=False, Run_Info=False)
            cReport.get_assay_info()
            cReport.get_testplate_info(WithStats=False,WithRunID=True)
            cReport.gen_pivot_tables()

            print(f" [Report] Project: {pk} [{cReport.n_compounds} {cReport.n_assays} {cReport.n_testplates} {cReport.n_screenruns} {cReport.n_sc} {cReport.n_dr}]")
            
            if cReport.n_samples>0:
                req = HttpResponse(content_type='application/vnd.ms-excel')
                req['Content-Disposition'] = f'attachment; filename={_xls_name}'
                cReport.to_excel(req)
                #return(req)
            #else:
            #return redirect(reverse("project_detail",kwargs={'pk':pk}))

    context={}
    context["object"]=_object

    return req

# -----------------------------------------------------------------
