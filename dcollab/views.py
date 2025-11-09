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
from django.http import JsonResponse, QueryDict
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy, reverse
from django.utils.functional import SimpleLazyObject
from django.views.generic import UpdateView, View

from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_RemoveView, Filtered_ListView, Htmx_UpdateView
 
# from applib.django.base.filters import Filtered_ListView
# from applib.django.base.views import permission_not_granted, Htmx_UpdateView, Base_CreateView, Base_UpdateView,  Base_RemoveView, File_CreateView

# from adjcoadd.constants import *

from dcollab.models import Organisation, Collab_Group, Collab_User
from dcollab.forms import (Organisation_Filter, Organisation_CreateForm, Organisation_UpdateForm,
                        CollabGroup_Filter, CollabGroup_CreateForm, CollabGroup_UpdateForm,
                        CollabGroup_Form, CollabMembership_FormSet,
                        CollabUser_Filter,  CollabUser_CreateForm,  CollabUser_UpdateForm,
                        )
# from dscreen.models import Screen_Run
# from dplate.models import MasterPlate, TestPlate
# from applib.report.screen_data import Report_Screening

#=================================================================================================
# Organisation
#=================================================================================================
class Organisation_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Organisation  
    template_name = 'dcollab/organisation/organisation_list.html'
    filterset_class = Organisation_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Organisation'
    app_name = 'dcollab'
    ordering=['-acreated_at']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = 'coadd_base.html'
        return context

# -----------------------------------------------------------------
@login_required
def Organisation_CreateView(req):
    '''
    View to Create new Organisation foreignkey: Dictionary. 
    '''
    print('Organisation_CreateView')  
    kwargs={}
    kwargs['user']=req.user
    form=Organisation_CreateForm()
    if req.method=='POST':
        form=Organisation_CreateForm(req.POST) 
        if form.is_valid():
            print('Organisation_CreateView Valid')
            try:
                with transaction.atomic(using='dcollab'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Organisation','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dcollab/organisation/organisation_create.html', { 'form':form, }) 

# -----------------------------------------------------------------
class Organisation_UpdateView(Htmx_UpdateView):

    form_class = Organisation_UpdateForm
    template_name = "dcollab/organisation/organisation_update.html"
    template_htmx = "dcollab/organisation/organisation_update_htmx.html"
    model = Organisation

#=================================================================================================
# Collaborator Group
#=================================================================================================
class CollabGroup_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Collab_Group  
    template_name = 'dcollab/collabgroup/collabgroup_list.html'
    filterset_class = CollabGroup_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Collab_Group'
    app_name = 'dcollab'
    ordering=['-acreated_at']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = 'coadd_base.html'
        return context

# -----------------------------------------------------------------
@login_required
def CollabGroup_CreateView(req):
    '''
    View to Create new Organisation foreignkey: Dictionary. 
    '''
    print('CollabGroup_CreateView')  
    kwargs={}
    kwargs['user']=req.user
    form=CollabGroup_CreateForm()
    if req.method=='POST':
        form=CollabGroup_CreateForm(req.POST) 
        if form.is_valid():
            print('CollabGroup_CreateView Valid')
            try:
                with transaction.atomic(using='dcollab'):
                    instance=form.save(commit=False) 
                    instance.save(**kwargs)
                    ApplicationLog.add('Create',str(instance.pk),'Info',req.user,str(instance.pk),'Create a new Collab Group','Completed')
                    return redirect(req.META['HTTP_REFERER'])
            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            messages.warning(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])          
    return render(req, 'dcollab/collabgroup/collabgroup_create.html', { 'form':form, }) 

# -----------------------------------------------------------------
@login_required
def CollabGroup_XUpdateView_old(req, pk):

    _object=get_object_or_404(Collab_Group, group_id=pk)

    kwargs={}
    kwargs['user']=req.user
    message={'status':'update','text':''}
    
    form=CollabGroup_UpdateForm(instance=_object)
    update_url = 'collabgroup_update'
    
    if req.method=='POST':
        try:
            with transaction.atomic(using='dscreen'):
                obj = Collab_Group.objects.select_for_update().get(group_id=pk)
                form= CollabGroup_UpdateForm(req.POST, instance=obj)    
                if form.is_valid():
                    instance=form.save(commit=False)
                    #update_screenrun_summary(instance)
                    instance.save(**kwargs)

                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update CollabGroup','Completed')
                    message={'status':'saved','text':f'CollabGroup [{pk}] Updated'}
                    return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':update_url, 'update_pk':pk}) 
                else:
                    messages.warning(req, f'Update failed due to {form.errors} error')
                    
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            message={'status':'update','text':'Input Error'}
            return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':update_url, 'update_pk':pk}) 

    context={}
    context["object"]=_object
    context["form"]=form
    
    return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':update_url, 'update_pk':pk}) 


#=================================================================================================
class CollabGroup_UpdateView(LoginRequiredMixin, UpdateView):
    
    model = Collab_Group
    form_class = CollabGroup_Form
    template_name = "dcollab/collab/group_update.html"
    success_url = reverse_lazy("collabgroup_list")

    # redirect unauthenticated users
    login_url = "login"              # your login URL name (can be path like '/accounts/login/')
    redirect_field_name = "next"     # optional, controls ?next= param

    def get_user_queryset(self):
        """Preload all users efficiently once."""
        return Collab_User.objects.select_related("organisation_id").order_by("last_name", "first_name")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        user_query_set = self.get_user_queryset()
        
        if self.request.POST:
            context["membership_formset"] = CollabMembership_FormSet(
                self.request.POST, instance=self.object
            )
        else:
            context["membership_formset"] = CollabMembership_FormSet(instance=self.object)
        return context

    def form_valid(self, form):
        context = self.get_context_data()
        membership_formset = context["membership_formset"]

        if membership_formset.is_valid():
            self.object = form.save()
            membership_formset.instance = self.object
            membership_formset.save()
            messages.success(self.request, "Group and members updated successfully.")
            if self.request.htmx:
                context["membership_formset"] = CollabMembership_FormSet(instance=self.object)
                return render(self.request, "dcollab/collab/partials/membership_list.html", context)
            return redirect(self.success_url)
        else:
            if self.request.htmx:
                return render(self.request, "dcollab/collab/partials/membership_list.html", context)
            return self.form_invalid(form)  
        
class AddMembershipRow_View(View):
    def get(self, request, pk):
        group = get_object_or_404(Collab_Group, pk=pk)
        formset = CollabMembership_FormSet(instance=group)
        new_form = formset.empty_form
        context = {"form": new_form}
        return render(request, "dcollab/collab/partials/membership_row.html", context)
             
#=================================================================================================
# Collaborator User
#=================================================================================================
class CollabUser_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Collab_User 
    template_name = 'dcollab/collabuser/collabuser_list.html'
    filterset_class = CollabUser_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Collab_User'
    app_name = 'dcollab'
    ordering=['-acreated_at']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['base_template'] = 'coadd_base.html'
        return context

# -----------------------------------------------------------------
@login_required
def CollabUser_CreateView(req):
    '''
    View to Create new ScreenRun foreignkey: Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    message={'status':'new','text':''}

    form=CollabUser_CreateForm()
    create_url = 'collabuser_create'
    
    #print(f" [CollabUser_CreateView] {req.method} {req.POST}")
    if req.method=='POST':
        form=CollabUser_CreateForm(req.POST) 
        if form.is_valid():
            #print(f" [CollabUser_CreateView] Valid Form")
            try:
                with transaction.atomic(using='dscreen'):
                    instance=form.save(commit=False)
                    instance.save(**kwargs) 
                    _newid = str(instance.user_id)
                    print(f" [CollabUser_CreateView] Saved:  [{_newid}]")            
                    message={'status':'saved','text':f'CollabUser [{_newid}] Created'}
                    return render(req, 'modal/createModel_partial_modal.html', {'message':message})

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    message={'status':'error','text':f'IntegrityError [{err}]'}
                    return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':create_url})
                    #return redirect(req.META['HTTP_REFERER'])                 
        else:
            messages.warning(req, form.errors)
            message={'status':'new','text':'Input Error'}
            return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':create_url})
            #return redirect(req.META['HTTP_REFERER'])          

    return render(req, 'modal/createModel_partial_modal.html', {'form':form, 'message':message, 'create_url':create_url}) 

# -----------------------------------------------------------------
@login_required
def CollabUser_UpdateView(req, pk):

    _object=get_object_or_404(Collab_User, assay_id=pk)

    kwargs={}
    kwargs['user']=req.user
    message={'status':'update','text':''}
    
    form=CollabUser_UpdateForm(instance=_object)
    
    if req.method=='POST':
        try:
            with transaction.atomic(using='dscreen'):
                obj = Collab_User.objects.select_for_update().get(assay_id=pk)
                form= CollabUser_UpdateForm(req.POST, instance=obj)    
                if form.is_valid():
                    instance=form.save(commit=False)
                    #update_screenrun_summary(instance)
                    instance.save(**kwargs)

                    ApplicationLog.add('Update',str(instance.pk),'Info',req.user,str(instance.pk),'Update CollabUser','Completed')
                    message={'status':'saved','text':f'CollabUser [{pk}] Updated'}
                    return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':'assay_update', 'update_pk':pk}) 
                else:
                    messages.warning(req, f'Update failed due to {form.errors} error')
                    
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            message={'status':'update','text':'Input Error'}
            return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':'assay_update', 'update_pk':pk}) 

    context={}
    context["object"]=_object
    context["form"]=form
    
    return render(req, 'modal/updateModel_partial_modal.html', {'form':form, 'message':message, 'update_url':'assay_update', 'update_pk':pk}) 
