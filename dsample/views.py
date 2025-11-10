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
from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_RemoveView, Filtered_ListView
 
# from applib.django.base.filters import Filtered_ListView
# from applib.django.base.views import permission_not_granted, Htmx_UpdateView, Base_CreateView, Base_UpdateView,  Base_RemoveView, File_CreateView

# from adjcoadd.constants import *

from dsample.models import Project, COADD_Compound, ABase_Compound_Batch
from dsample.forms import Project_Filter, Project_CreateForm, Project_UpdateForm
from dsample.utils.summary import update_project_summary
from dscreen.models import Screen_Run
from dplate.models import MasterPlate, TestPlate
from applib.report.screen_data import Report_Screening
from applib.project.stockprep_project import StockPrep_Project


#=================================================================================================
# Project
#=================================================================================================
class Project_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Project  
    template_name = 'dsample/project/project_list.html'
    filterset_class = Project_Filter
    model_fields = model.LIST_VIEW_FIELDS
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
    View to Create new Project  
    '''  
    
    kwargs={}
    kwargs['user']=req.user
    message={'status':'new','text':''}

    form=Project_CreateForm()
    
    #print(f" [Project_CreateView] {req.method} {req.POST}")
    if req.method=='POST':
        form=Project_CreateForm(req.POST) 
        if form.is_valid():
            #print(f" [Project_CreateView] Valid Form")
            try:
                with transaction.atomic(using='dsample'):
                    instance=form.save(commit=False)
                    instance.save(**kwargs) 
                    _newid = str(instance.project_id)
                    print(f" [Project_CreateView] Saved:  [{_newid}]")            
                    message={'status':'saved','text':f'Project [{_newid}] Created'}
                    return render(req, 'modal/createModel_partial_modal.html', {'message':message})

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    message={'status':'error','text':f'IntegrityError [{err}]'}
                    return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':'project_create'})
                    #return redirect(req.META['HTTP_REFERER'])                 
        else:
            messages.warning(req, form.errors)
            message={'status':'new','text':'Input Error'}
            return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':'project_create'})
            #return redirect(req.META['HTTP_REFERER'])          

    return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':'project_create'}) 

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


    # Paginated and filtered list
    paginate_by = 50
    _coadd_compounds = COADD_Compound.objects.filter(project_id=_object, )
    context["coadd_objs"] = _coadd_compounds
    context["coadd_count"] = _coadd_compounds.count()
    context["coadd_fields"] = COADD_Compound.get_fields()

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
                    update_project_summary(instance)
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
class Project_RemoveView(Base_RemoveView):
    model = Project
    transaction_use = 'dsample'

# -----------------------------------------------------------------
@login_required
def Project_ReportView(req, pk):
    print(f" [Report] Project: {req}")
    
    if req.method=='GET':

        _object=get_object_or_404(Project, project_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'Project_{pk}_Summary_{_now:%Y%m%d}.xlsx'

        cReport = Report_Screening()
        cReport.qry_by_ProjectID(_object)
        if cReport.n_compounds>0:
            cReport.get_dataframe()
            cReport.get_sample_info(Storage_Info=False, Structure_Info=True, Run_Info=False)
            cReport.get_assay_info()
            cReport.get_testplate_info(WithStats=False,WithRunID=True)
            cReport.gen_pivot_tables(PivRows=['sample_class','sample_code','sample_id'])

            print(f" [Report] Project: {pk} [{cReport.n_compounds} {cReport.n_assays} {cReport.n_testplates} {cReport.n_screenruns} {cReport.n_sc} {cReport.n_dr}]")
            
            if cReport.n_samples>0:
                print(f" [Report] Project: {cReport.n_samples}")
                req = HttpResponse(content_type='application/vnd.ms-excel')
                req['Content-Disposition'] = f'attachment; filename={_xls_name}'
                print(f" [Report] Project: {req['Content-Disposition']}")
                cReport.to_excel(req)
                print(f" [Report] Project: done {req}")

        return(req)
    
# -----------------------------------------------------------------

@login_required
def Project_StockPrepView(req, pk):

    if req.method=='GET':

        _object=get_object_or_404(Project, project_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'Project_{pk}_CmpdPrep_{_now:%Y%m%d}.xlsx'

        cCmpdPrep = StockPrep_Project(pk)
        cCmpdPrep.get_samples()

        if cCmpdPrep.n_samples>0:
            req = HttpResponse(content_type='application/vnd.ms-excel')
            req['Content-Disposition'] = f'attachment; filename={_xls_name}'
            cCmpdPrep.to_excel(req)
        return req
