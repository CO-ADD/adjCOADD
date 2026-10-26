import os
import json
import datetime
from io import BytesIO as IO

from rdkit import Chem
from django_filters.views import FilterView

from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.conf import settings
from django.contrib import messages
from django.core.exceptions import FieldDoesNotExist, ValidationError
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

from dsample.models import Project, COADD_Compound, ABase_Compound_Batch, Project_Membership, ABase_Compound
from dsample.forms import Project_Filter, Project_CreateForm, Project_CreateMinimalForm, Project_UpdateForm, ABase_Filter, COADD_Filter
from dsample.utils.summary import update_project_summary
from dscreen.models import Screen_Run
from dplate.models import MasterPlate, TestPlate
from applib.report.screen_data import Report_Screening
from applib.project.stockprep_project import StockPrep_Project
from dcollab.models import Collab_User, Collab_Group


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

    _empty_group = Collab_Group.get(Project.EMPTY_GROUP)
    form=Project_CreateMinimalForm()
    
    
    #print(f" [Project_CreateView] {req.method} {req.POST}")
    if req.method=='POST':
        form=Project_CreateMinimalForm(req.POST) 
        if form.is_valid():
            #print(f" [Project_CreateView] Valid Form")
            try:
                with transaction.atomic(using='dsample'):
                    print(f" [Project_CreateView] Saving Project")   
                    
                    instance=form.save(commit=False)
                    
                    instance.set_defaults_model()
                    instance.group_id = _empty_group
                    
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
def Project_CreateMinimalView(req):
    '''
    View to Create new Project  
    '''  
    
    kwargs={}
    kwargs['user']=req.user
    message={'status':'new','text':''}

    _empty_user = Collab_User.get(Project.EMPTY_USER)
    _empty_group = Collab_Group.get(Project.EMPTY_GROUP)
    
    #form=Project_CreateMinimalForm(initial={'group_id':_empty_group, 'project_members':[_empty_user] })
    #form=Project_CreateMinimalForm(initial={'group_id':_empty_group, })

    form=Project_CreateMinimalForm()
    if req.method=='POST':
        form=Project_CreateMinimalForm(req.POST) 
        if form.is_valid():
            try:
                with transaction.atomic(using='dsample'):

                    print(f" [Project_CreateView] Saving Project")   
                    instance=form.save(commit=False)
                    
                    print(f" [Project_CreateView] Set Defaults {instance.project_id} {instance.project_name} {instance.project_type} {instance.project_status} {instance.pub_status} {instance.provided_container}")  
                    #print(instance.project_members.all())

                    #instance.set_defaults_model()
                    #instance.set_none_field('ora_project_id')
                                
                    #print(f" [Project_CreateView] Save Project {kwargs}")
                    #valDict = instance.validate_model()
                    #print(instance.project_members)
                    #print(f" [Project_CreateView] Validation {valDict}")
                    
                    instance.save(verbose=1)
                    print(f" [Project_CreateView] Saved")
                    form.save_m2m() 

                    print(f" [Project_CreateView] Init Group")
                    #instance.init_group()

                    print(f" [Project_CreateView] Wrapup")
                    _newid = str(instance.project_id)
                    message={'status':'saved','text':f'Project [{_newid}] Created'}
                    print(f" [Project_CreateView] Saved:  [{_newid}]")
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
    _object=get_object_or_404(Project, project_id=pk)
    context["object"]=_object
    context["display_groups"] = _project_display_groups(_object)

    # Paginated and filtered list
    _view = req.GET.get('view', 'table')
    _per_page = 15 if _view == 'card' else 50

    _coadd_compounds = COADD_Compound.objects.filter(project_id=_object, astatus__gt=0)
    _coadd_paginator = Paginator(_coadd_compounds, _per_page)
    _coadd_page = _coadd_paginator.get_page(req.GET.get('page'))

    context["view"] = _view
    context["coadd_objs"] = _coadd_page
    context["coadd_page_obj"] = _coadd_page
    context["coadd_count"] = _coadd_compounds.count()
    context["coadd_fields"] = COADD_Compound.get_fields()

    _collab = _object.get_members()
    context["collab_objs"] = _collab
    context["collab_count"] = len(_collab)
    context["collab_fields"] = ['Name','Role','Status']

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
def _project_display_groups(instance):
    """
    Build read-only label/value rows for the Project detail view,
    grouped according to Project.DETAIL_VIEW_GROUPS.
    """
    groups = []
    for field_names in Project.DETAIL_VIEW_GROUPS:
        rows = []
        for name in field_names:
            try:
                field = Project._meta.get_field(name)
            except FieldDoesNotExist:
                continue
            value = getattr(instance, name, None)
            if value is None or value == '':
                value = ''
            elif isinstance(value, (list, tuple)):
                value = ', '.join(str(v) for v in value)
            else:
                value = str(value)
            rows.append({'label': field.verbose_name, 'value': value})
        groups.append(rows)
    return groups

# -----------------------------------------------------------------
@login_required
def Project_UpdateView(req, pk):
    _object=get_object_or_404(Project, project_id=pk)

    if req.method=='POST':
        form=Project_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dsample'):
                    instance=form.save(commit=False)
                    update_project_summary(instance)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update Project','Completed')
                return redirect('project_detail', pk=pk)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('project_detail', pk=pk)

    # GET - render the editable form for the modal
    form=Project_UpdateForm(instance=_object)
    context={}
    context["object"]=_object
    context["form"]=form
    return render(req, "dsample/project/project_update_modal.html", context)

# -----------------------------------------------------------------
@login_required
def Project_MembersEditView(req, pk):
    """
    Manage the project_members (PC/AC/VC). The Lead Investigator (LI) is
    excluded - it is derived from Project.group_id.li_user_id.
    """
    _object=get_object_or_404(Project, project_id=pk)

    if req.method == 'POST':
        _roles = dict(Project_Membership.MEMBERSHIP_CHOICES)

        # 1) Remove selected memberships
        for user_id in req.POST.getlist('remove'):
            Project_Membership.objects.filter(project_id=_object, user_id=user_id).delete()

        # 2) Change roles for remaining memberships
        for _memb in Project_Membership.objects.filter(project_id=_object):
            _role = req.POST.get(f'role_{_memb.user_id_id}')
            if _role in _roles and _role != _memb.role:
                _memb.role = _role
                _memb.save()

        # 3) Add a new member
        _add_user = req.POST.get('add_user')
        _add_role = req.POST.get('add_role')
        if _add_user and _add_role in _roles:
            _user = Collab_User.objects.filter(user_id=_add_user).first()
            if _user:
                _memb = Project_Membership.objects.filter(project_id=_object, user_id=_user).first()
                if _memb:
                    _memb.role = _add_role
                    _memb.save()
                else:
                    Project_Membership.objects.create(project_id=_object, user_id=_user, role=_add_role, status='C')

        ApplicationLog.add('Update', 'Project_Members', 'Info', req.user, str(_object.pk), 'Update Project Members', 'Completed')
        return redirect('project_detail', pk=pk)

    # GET - render the collaborator edit modal
    _memberships = Project_Membership.objects.filter(project_id=_object).select_related('user_id')
    _member_ids = list(_memberships.values_list('user_id_id', flat=True))
    _available = Collab_User.objects.exclude(user_id__in=_member_ids).order_by('first_name', 'last_name')

    context = {
        'object': _object,
        'memberships': _memberships,
        'roles': Project_Membership.MEMBERSHIP_CHOICES,
        'available_users': _available,
    }
    return render(req, "dsample/project/project_members_modal.html", context)

# -----------------------------------------------------------------
class Project_RemoveView(Base_RemoveView):
    model = Project
    transaction_use = 'dsample'

# -----------------------------------------------------------------
@login_required
def Project_ReportView(req, pk):
    _object=get_object_or_404(Project, project_id=pk)

    if req.method == 'POST':
        # Collect the selected report options from the modal form
        options = {
            'Storage_Info':   'Storage_Info' in req.POST,
            'Structure_Info': 'Structure_Info' in req.POST,
            'Run_Info':       'Run_Info' in req.POST,
            'WithStats':      'WithStats' in req.POST,
            'WithRunID':      'WithRunID' in req.POST,
        }
        report_name, report_url = _generate_project_report(pk, options)

        if report_url:
            message = {
                'status': 'generated',
                'report_name': report_name,
                'report_url': report_url,
            }
        else:
            message = {
                'status': 'error',
                'text': 'No screening data available for this project.',
            }
        return render(req, 'dsample/project/project_report_modal.html', {'message': message, 'object': _object})

    # GET - render the report options for the modal
    options = {
        'Storage_Info':   False,
        'Structure_Info': True,
        'Run_Info':       False,
        'WithStats':      False,
        'WithRunID':      True,
    }
    message = {'status': 'options'}
    return render(req, 'dsample/project/project_report_modal.html', {'message': message, 'options': options, 'object': _object})


def _generate_project_report(pk, options):
    """
    Generate the Project screening report xlsx using the selected options.
    Saves the report under MEDIA_ROOT/reports/ and returns (report_name, report_url).
    Returns (None, None) when there is no screening data to report.
    """
    _object = get_object_or_404(Project, project_id=pk)
    _now = datetime.datetime.now()
    report_name = f'Project_{pk}_Summary_{_now:%Y%m%d%H%M%S}.xlsx'

    cReport = Report_Screening()
    cReport.qry_by_ProjectID(_object)

    if cReport.n_compounds == 0:
        return None, None

    cReport.get_dataframe()
    cReport.get_sample_info(Storage_Info=options['Storage_Info'],
                            Structure_Info=options['Structure_Info'],
                            Run_Info=options['Run_Info'])
    cReport.get_assay_info()
    cReport.get_testplate_info(WithStats=options['WithStats'],
                               WithRunID=options['WithRunID'])
    cReport.gen_pivot_tables(PivRows=['sample_class','sample_code','sample_id'])

    print(f" [Report] Project: {pk} [{cReport.n_compounds} {cReport.n_assays} {cReport.n_testplates} {cReport.n_screenruns} {cReport.n_sc} {cReport.n_dr}]")

    if cReport.n_samples == 0:
        return None, None

    report_dir = os.path.join(settings.MEDIA_ROOT, 'reports')
    os.makedirs(report_dir, exist_ok=True)
    report_path = os.path.join(report_dir, report_name)
    cReport.to_excel(report_path)

    report_url = f"/{settings.MEDIA_URL.strip('/')}/reports/{report_name}"
    return report_name, report_url
    
# -----------------------------------------------------------------

@login_required
def Project_StockPrepView(req, pk):

    if req.method=='GET':

        _object=get_object_or_404(Project, project_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'{pk}_StockPrep_{_now:%d%b%Y}.xlsx'

        cCmpdPrep = StockPrep_Project(pk)
        cCmpdPrep.get_samples()

        if cCmpdPrep.n_samples>0:
            cCmpdPrep.make_stockprep()
            req = HttpResponse(content_type='application/vnd.ms-excel')
            req['Content-Disposition'] = f'attachment; filename={_xls_name}'
            cCmpdPrep.to_excel(req)
        return req

#=================================================================================================
# ABase Compound
#=================================================================================================

ABase_DETAIL_FIELDS = ['compound_id','compound_code','compound_name','compound_desc',
                       'reg_mw','reg_mf','structure_type','structure_metal','structure_id']


class ABase_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = ABase_Compound
    template_name = 'dsample/abase/abase_list.html'
    filterset_class = ABase_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'ABase_Compound'
    app_name = 'dsample'
    ordering = ['compound_id']
    paginate_by = 50

    def get_queryset(self):
        # Filter-first: do not evaluate the (potentially huge) queryset until
        # at least one filter/search term has been supplied.
        if not self.get_filter_record():
            self.filterset = self.filterset_class(self.request.GET,
                                                  queryset=self.model.objects.none(),
                                                  filterset_dict={}, deep=False)
            self.filter_count = 0
            return self.model.objects.none()
        return super().get_queryset()

    def get_context_data(self, **kwargs):
        # Replicate Filtered_ListView.get_context_data() but avoid the
        # full-table COUNT() (expensive for very large tables).
        context = super(Filtered_ListView, self).get_context_data(**kwargs)

        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        context['filter'] = self.filterset
        context['paginate_by'] = self.request.GET.get("paginate_by", self.paginate_by)
        context['fields'] = self.model.get_fields(fields=self.model_fields)
        context['filterset'] = filter_record
        context['querycount'] = self.filter_count
        context['Count'] = None
        context['has_filters'] = bool(filter_record_dict)
        return context


# -----------------------------------------------------------------
@login_required
def ABase_DetailView(req, pk):
    """
    Detail view for an ABase_Compound entry:
      - read-only field groups
      - chemical structure rendering (reg_smiles) on the left
      - subsection accordion with the related ABase_Compound_Batch entries
    """
    _object = get_object_or_404(ABase_Compound, compound_id=pk)

    _rows = []
    for name in ABase_DETAIL_FIELDS:
        try:
            field = ABase_Compound._meta.get_field(name)
        except FieldDoesNotExist:
            continue
        value = getattr(_object, name, None)
        if value is None or value == '':
            value = ''
        else:
            value = str(value)
        _rows.append({'label': field.verbose_name, 'value': value})

    # Split rows into column groups for the display template
    _groups = []
    _chunk = 7
    for i in range(0, len(_rows), _chunk):
        _groups.append(_rows[i:i+_chunk])

    _batches = ABase_Compound_Batch.objects.filter(compound_id=_object)

    context = {
        'object': _object,
        'display_groups': _groups,
        'batch_objs': _batches,
        'batch_count': _batches.count(),
        'batch_fields': ABase_Compound_Batch.get_fields(),
    }
    return render(req, "dsample/abase/abase_detail.html", context)


#=================================================================================================
# COADD Compound
#=================================================================================================

COADD_DETAIL_FIELDS = ['compound_id','compound_code','compound_name','compound_desc',
                       'project_id','compound_type','compound_source',
                       'std_status','std_process','std_issues','std_nfrag',
                       'std_salt','std_ion','std_solvent','std_metal','std_structure_type',
                       'std_mw','std_mf','pub_status','pub_date']


class COADD_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = COADD_Compound
    template_name = 'dsample/coadd/coadd_list.html'
    filterset_class = COADD_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'COADD_Compound'
    app_name = 'dsample'
    ordering = ['compound_id']
    paginate_by = 50

    def get_queryset(self):
        # Filter-first: do not evaluate the (very large) queryset until
        # at least one filter/search term has been supplied.
        if not self.get_filter_record():
            self.filterset = self.filterset_class(self.request.GET,
                                                  queryset=self.model.objects.none(),
                                                  filterset_dict={}, deep=False)
            self.filter_count = 0
            return self.model.objects.none()
        return super().get_queryset()

    def get_context_data(self, **kwargs):
        # Replicate Filtered_ListView.get_context_data() but avoid the
        # full-table COUNT() (expensive for very large tables).
        context = super(Filtered_ListView, self).get_context_data(**kwargs)

        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        context['filter'] = self.filterset
        context['paginate_by'] = self.request.GET.get("paginate_by", self.paginate_by)
        context['fields'] = self.model.get_fields(fields=self.model_fields)
        context['filterset'] = filter_record
        context['querycount'] = self.filter_count
        context['Count'] = None
        context['has_filters'] = bool(filter_record_dict)
        return context


# -----------------------------------------------------------------
@login_required
def COADD_DetailView(req, pk):
    """
    Detail view for a COADD_Compound entry:
      - read-only field groups (no subtables for now)
      - chemical structure rendering from std_smiles on the left
    """
    _object = get_object_or_404(COADD_Compound, compound_id=pk)

    _rows = []
    for name in COADD_DETAIL_FIELDS:
        try:
            field = COADD_Compound._meta.get_field(name)
        except FieldDoesNotExist:
            continue
        value = getattr(_object, name, None)
        if value is None or value == '':
            value = ''
        else:
            value = str(value)
        _rows.append({'label': field.verbose_name, 'value': value})

    # Split rows into column groups for the display template
    _groups = []
    _chunk = 7
    for i in range(0, len(_rows), _chunk):
        _groups.append(_rows[i:i+_chunk])

    context = {
        'object': _object,
        'display_groups': _groups,
    }
    return render(req, "dsample/coadd/coadd_detail.html", context)


