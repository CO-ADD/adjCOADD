
import os
import json
import datetime

#from rdkit import Chem
#from django_filters.views import FilterView

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
from django.utils.safestring import mark_safe

from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_RemoveView, Filtered_ListView

# from adjcoadd.constants import *

from dscreen.models import Screen_Run
#from dsummary.models import Summary_ScreenRun 
from dscreen.forms import ScreenRun_Filter, ScreenRun_CreateForm, ScreenRun_UpdateForm
from dscreen.utils.screenrun_process import Upload_ReadOuts_Process
from dscreen.utils.summary import update_screenrun_summary, get_projects_screenrun
from dsample.models import Project
from dplate.models import MasterPlate, TestPlate, TestWell
from applib.report.screen_data import Report_Screening

#=================================================================================================
# ScreenRun
#=================================================================================================
class ScreenRun_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Screen_Run  
    template_name = 'dscreen/screenrun/screenrun_list.html'
    filterset_class = ScreenRun_Filter
    model_fields = model.LIST_VIEW_FIELDS
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
    # if req.method == 'GET':
    #     print(f"[ScreenRun_DetailView] GET {req.GET}")
    # if req.method == 'POST':
    #     print(f"[ScreenRun_DetailView] POST: {req.POST}")

    context["object"]=_object
    context["form"]=form

    if str(_object.run_type) in ['HCR','PSR']:
        context["process"] = {"type":"Screening"}
        
        # Paginated and filtered list
        paginate_by = 50
        _testplates = TestPlate.objects.filter(run_id=_object, )
        context["testplate_objs"] = _testplates
        context["testplate_count"] = _testplates.count()
        context["testplate_fields"]=TestPlate.get_fields()
        
        # Single page list
        _prj_dict = get_projects_screenrun(_object)
        _projects = [p['project_id'] for p in _prj_dict]
        _projects = Project.objects.filter(project_id__in=_projects, )
        context["project_objs"] = _projects
        context["project_count"] = _projects.count()
        context["project_fields"]= Project.get_fields()
        
    elif str(_object.run_type) in ['SEQ']:
        context["process"] = {"type":"Sequencing"}
    else:
        context["process"] = {"type":"Undefined"}
        
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
    if req.method=='POST':
        try:
            with transaction.atomic(using='dscreen'):
                obj = Screen_Run.objects.select_for_update().get(run_id=pk)
                form= ScreenRun_UpdateForm(req.POST, instance=obj)    
                if form.is_valid():
                    instance=form.save(commit=False)
                    update_screenrun_summary(instance)
                    instance.save(**kwargs)

                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update Screen_Run','Completed')

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
class ScreenRun_RemoveView(Base_RemoveView):
    model = Screen_Run
    transaction_use = 'dscreen'


# -----------------------------------------------------------------
@login_required
def ScreenRun_ReportView(req, pk):

    if req.method=='GET':

        _object=get_object_or_404(Screen_Run, run_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'ScreenRun_{pk}_Summary_{_now:%Y%m%d}.xlsx'

        cReport = Report_Screening()
        cReport.qry_by_RunID([pk])
        cReport.get_dataframe()
        cReport.get_sample_info(Storage_Info=True,Structure_Info=True)
        cReport.get_assay_info()
        cReport.get_testplate_info()
        cReport.add_hcr_selection()

        #cReport.add_vitek_ast()
        #cReport.add_antibiogram_data(cReport.ORGANISMS['COADD'])

        cReport.gen_pivot_tables(PivTables = ['Values','AssayID','Act'], PivRows=['project_id','sample_class','sample_code','sample_id'])

        if cReport.n_samples>0:
            req = HttpResponse(content_type='application/vnd.ms-excel')
            req['Content-Disposition'] = f'attachment; filename={_xls_name}'
            cReport.to_excel(req)
        return req
