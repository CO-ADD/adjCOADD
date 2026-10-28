import os
import json
from rdkit import Chem

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import FieldDoesNotExist, ValidationError
from django.core.paginator import Paginator
from django.db import transaction, IntegrityError
from django.shortcuts import get_object_or_404, render, redirect

from applib.django.base.views import Base_RemoveView, Filtered_ListView
from apputil.models import ApplicationLog
from apputil.forms import Document_Form

from dpeptide.models import (Peptide, Peptide_Isoform, Isoform_Batch,
                             AminoAcid_Dictionary, Modification_Dictionary)
from dpeptide.forms import (Peptide_Filter, Peptide_CreateForm, Peptide_UpdateForm,
                            Isoform_Filter, Isoform_CreateForm, Isoform_UpdateForm,
                            IsoBatch_CreateForm, IsoBatch_UpdateForm,
                            AminoAcid_Filter, AminoAcid_Form, AminoAcid_UpdateForm,
                            Modification_Filter, Modification_Form, Modification_UpdateForm)

#=================================================================================================
# Helpers
#=================================================================================================
def _display_groups(instance, group_defs):
    """Build read-only label/value row groups for a detail view."""
    groups = []
    for field_names in group_defs:
        rows = []
        for name in field_names:
            try:
                field = instance._meta.get_field(name)
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
    ordering = ['-acreated_at']


# -----------------------------------------------------------------
@login_required
def Peptide_CreateView(req):
    form = Peptide_CreateForm()

    if req.method == 'POST':
        form = Peptide_CreateForm(req.POST)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create a new Peptide', 'Completed')
                return redirect('peptide_list')
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('peptide_list')

    context = {'form': form, 'modal_title': 'Create Peptide', 'form_action': 'peptide_create'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
@login_required
def Peptide_DetailView(req, pk):
    _object = get_object_or_404(Peptide, peptide_id=pk)
    context = {
        'object': _object,
        'display_groups': _display_groups(_object, Peptide.DETAIL_VIEW_GROUPS),
        'doc_form': Document_Form,
    }
    return render(req, 'dpeptide/peptide/peptide_detail.html', context)


# -----------------------------------------------------------------
@login_required
def Peptide_UpdateView(req, pk):
    _object = get_object_or_404(Peptide, peptide_id=pk)

    if req.method == 'POST':
        form = Peptide_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Peptide', 'Completed')
                return redirect('peptide_detail', pk=pk)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('peptide_detail', pk=pk)

    form = Peptide_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Peptide', 'modal_subtitle': pk, 'form_action': 'peptide_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
class Peptide_RemoveView(Base_RemoveView):
    model = Peptide
    transaction_use = 'dpeptide'


#=================================================================================================
# Peptide Isoform
#=================================================================================================
class Isoform_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Peptide_Isoform
    template_name = 'dpeptide/isoform/isoform_list.html'
    filterset_class = Isoform_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Peptide_Isoform'
    app_name = 'dpeptide'
    ordering = ['isoform_id']


# -----------------------------------------------------------------
@login_required
def Isoform_CreateView(req):
    form = Isoform_CreateForm(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create a new Isoform', 'Completed')
                return redirect('isoform_detail', pk=instance.pk)
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('isoform_list')

    context = {'form': form, 'modal_title': 'Create Isoform', 'form_action': 'isoform_create'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
@login_required
def Isoform_DetailView(req, pk):
    _object = get_object_or_404(Peptide_Isoform, isoform_id=pk)
    _peptide = _object.peptide_id

    context = {
        'object': _object,
        'peptide': _peptide,
        'display_groups': _display_groups(_object, Peptide_Isoform.DETAIL_VIEW_GROUPS),
        'peptide_groups': _display_groups(_peptide, Peptide.DETAIL_VIEW_GROUPS),
        'doc_form': Document_Form,
    }

    # Isoform batches - paginated
    _batches = Isoform_Batch.objects.filter(isoform_id=_object, astatus__gte=0)
    _paginator = Paginator(_batches, 25)
    _batch_page = _paginator.get_page(req.GET.get('page'))

    context['batch_objs'] = _batch_page
    context['batch_page_obj'] = _batch_page
    context['batch_count'] = _batches.count()
    context['batch_fields'] = Isoform_Batch.get_fields()

    return render(req, 'dpeptide/isoform/isoform_detail.html', context)


# -----------------------------------------------------------------
@login_required
def Isoform_UpdateView(req, pk):
    _object = get_object_or_404(Peptide_Isoform, isoform_id=pk)

    if req.method == 'POST':
        form = Isoform_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Isoform', 'Completed')
                return redirect('isoform_detail', pk=pk)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('isoform_detail', pk=pk)

    form = Isoform_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Isoform', 'modal_subtitle': pk, 'form_action': 'isoform_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
class Isoform_RemoveView(Base_RemoveView):
    model = Peptide_Isoform
    transaction_use = 'dpeptide'


#=================================================================================================
# Isoform Batch
#=================================================================================================
@login_required
def IsoBatch_CreateView(req, isoform_id):
    _isoform = get_object_or_404(Peptide_Isoform, isoform_id=isoform_id)
    form = IsoBatch_CreateForm(req.POST or None)

    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.isoform_id = _isoform
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create a new Isoform Batch', 'Completed')
                return redirect('isoform_detail', pk=isoform_id)
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('isoform_detail', pk=isoform_id)

    context = {
        'form': form, 'object': _isoform,
        'modal_title': 'New Isoform Batch', 'modal_subtitle': isoform_id,
        'form_action': 'isobatch_create',
    }
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
@login_required
def IsoBatch_UpdateView(req, pk):
    _object = get_object_or_404(Isoform_Batch, isoformbatch_id=pk)

    if req.method == 'POST':
        form = IsoBatch_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Isoform Batch', 'Completed')
                return redirect('isoform_detail', pk=_object.isoform_id.isoform_id)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('isoform_detail', pk=_object.isoform_id.isoform_id)

    form = IsoBatch_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Isoform Batch', 'modal_subtitle': pk, 'form_action': 'isobatch_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


# -----------------------------------------------------------------
class IsoBatch_RemoveView(Base_RemoveView):
    model = Isoform_Batch
    transaction_use = 'dpeptide'


#=================================================================================================
# AminoAcid Dictionary
#=================================================================================================
class AminoAcid_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = AminoAcid_Dictionary
    template_name = 'dpeptide/dict/dictionary_list.html'
    filterset_class = AminoAcid_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'AminoAcid_Dictionary'
    app_name = 'dpeptide'
    ordering = ['aa_code']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Amino Acids'
        context['model_name'] = self.model_name
        context['list_url'] = 'aminoacid_list'
        context['create_url'] = 'aminoacid_create'
        context['update_url'] = 'aminoacid_update'
        return context


@login_required
def AminoAcid_CreateView(req):
    form = AminoAcid_Form(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create Amino Acid', 'Completed')
                return redirect('aminoacid_list')
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('aminoacid_list')

    context = {'form': form, 'modal_title': 'Create Amino Acid', 'form_action': 'aminoacid_create'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


@login_required
def AminoAcid_UpdateView(req, pk):
    _object = get_object_or_404(AminoAcid_Dictionary, aa_id=pk)
    if req.method == 'POST':
        form = AminoAcid_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Amino Acid', 'Completed')
                return redirect('aminoacid_list')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('aminoacid_list')

    form = AminoAcid_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Amino Acid', 'modal_subtitle': pk, 'form_action': 'aminoacid_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


#=================================================================================================
# Modification Dictionary
#=================================================================================================
class Modification_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Modification_Dictionary
    template_name = 'dpeptide/dict/dictionary_list.html'
    filterset_class = Modification_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Modification_Dictionary'
    app_name = 'dpeptide'
    ordering = ['pepmod_code']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Modifications'
        context['model_name'] = self.model_name
        context['list_url'] = 'modification_list'
        context['create_url'] = 'modification_create'
        context['update_url'] = 'modification_update'
        return context


@login_required
def Modification_CreateView(req):
    form = Modification_Form(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create Modification', 'Completed')
                return redirect('modification_list')
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('modification_list')

    context = {'form': form, 'modal_title': 'Create Modification', 'form_action': 'modification_create'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


@login_required
def Modification_UpdateView(req, pk):
    _object = get_object_or_404(Modification_Dictionary, pepmod_id=pk)
    if req.method == 'POST':
        form = Modification_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Modification', 'Completed')
                return redirect('modification_list')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('modification_list')

    form = Modification_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Modification', 'modal_subtitle': pk, 'form_action': 'modification_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)
