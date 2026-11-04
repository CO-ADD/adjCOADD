import os
import json
from rdkit import Chem

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import FieldDoesNotExist, ValidationError
from django.core.paginator import Paginator
from django.db import transaction, IntegrityError
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render, redirect

from applib.django.base.views import Base_RemoveView, Filtered_ListView
from applib.bio.biln import render_biln_svg, render_helm_svg
from apputil.models import ApplicationLog
from apputil.forms import Document_Form

from dpeptide.models import Peptide, Peptide_Isoform, Isoform_Batch, Modification_Dictionary, Isoform_Modification, Monomer
from dpeptide.forms import (Peptide_Filter, Peptide_CreateForm, Peptide_UpdateForm,
                            Isoform_Filter, Isoform_CreateForm, Isoform_UpdateForm,
                            IsoBatch_CreateForm, IsoBatch_UpdateForm,
                            Modification_Filter, Modification_Form, Modification_UpdateForm,
                            IsoformModification_Form,
                            Monomer_Filter, Monomer_Form, Monomer_UpdateForm)

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


def _error_svg(message):
    """Minimal SVG showing an error message, for <img> previews."""
    text = message.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="500" height="80">'
        f'<rect width="100%" height="100%" fill="#fff3cd" stroke="#ffc107"/>'
        f'<text x="10" y="25" font-family="sans-serif" font-size="12" fill="#664d03">{text}</text>'
        f'</svg>'
    )


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

    # Isoform modifications - grouped by position, ordered by mod_order
    _mods = _object.modification_links.select_related('pepmod_id').order_by('mod_position', 'mod_order')
    _mod_groups = []
    for pos_value, pos_label in Isoform_Modification.MOD_POSITION_CHOICES:
        _group_links = [link for link in _mods if link.mod_position == pos_value]
        if _group_links:
            _mod_groups.append({'position': pos_label, 'links': _group_links})
    context['modification_groups'] = _mod_groups
    context['modification_count'] = _mods.count()

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
# Isoform Modifications (through-model links)
#=================================================================================================

@login_required
def IsoformModification_CreateView(req, pk):
    _isoform = get_object_or_404(Peptide_Isoform, isoform_id=pk)
    form = IsoformModification_Form(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    link = form.save(commit=False)
                    link.isoform_id = _isoform
                    link.save()
                    ApplicationLog.add('Create', f"{_isoform.pk}:{link.pepmod_id}",
                                       'Info', req.user, str(_isoform.pk),
                                       f'Add modification {link.pepmod_id} ({link.mod_position}:{link.mod_order})', 'Completed')
                return redirect('isoform_detail', pk=pk)
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} - modification with same position/order may already exist.')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('isoform_detail', pk=pk)

    context = {'form': form, 'modal_title': 'Add Modification', 'modal_subtitle': pk}
    return render(req, 'dpeptide/isoform/isoform_modification_form_modal.html', context)


@login_required
def IsoformModification_UpdateView(req, pk):
    _link = get_object_or_404(Isoform_Modification, pk=pk)
    if req.method == 'POST':
        form = IsoformModification_Form(req.POST, instance=_link)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    form.save()
                    ApplicationLog.add('Update', f"{_link.isoform_id_id}:{_link.pepmod_id}",
                                       'Info', req.user, str(_link.isoform_id_id),
                                       'Update isoform modification', 'Completed')
                return redirect('isoform_detail', pk=_link.isoform_id_id)
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} - modification with same position/order may already exist.')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('isoform_detail', pk=_link.isoform_id_id)

    form = IsoformModification_Form(instance=_link)
    context = {'form': form, 'modal_title': 'Update Modification', 'modal_subtitle': _link.isoform_id_id}
    return render(req, 'dpeptide/isoform/isoform_modification_form_modal.html', context)


@login_required
def IsoformModification_RemoveView(req, pk):
    _link = get_object_or_404(Isoform_Modification, pk=pk)
    _isoform_pk = _link.isoform_id_id
    if req.method == 'POST':
        with transaction.atomic(using='dpeptide'):
            ApplicationLog.add('Remove', f"{_isoform_pk}:{_link.pepmod_id}",
                               'Warning', req.user, str(_isoform_pk),
                               f'Remove modification {_link.pepmod_id} ({_link.mod_position}:{_link.mod_order})', 'Completed')
            _link.delete()
        messages.success(req, f'Removed modification {_link.pepmod_id}')
    return redirect('isoform_detail', pk=_isoform_pk)


@login_required
def IsoformModification_ReorderView(req, pk, direction):
    _link = get_object_or_404(Isoform_Modification, pk=pk)
    _isoform_pk = _link.isoform_id_id
    if req.method == 'POST':
        _siblings = Isoform_Modification.objects.filter(
            isoform_id=_link.isoform_id, mod_position=_link.mod_position).order_by('mod_order')
        if direction == 'up':
            _other = _siblings.filter(mod_order__lt=_link.mod_order).order_by('-mod_order').first()
        else:
            _other = _siblings.filter(mod_order__gt=_link.mod_order).order_by('mod_order').first()

        if _other is not None:
            # mod_order is PositiveSmallIntegerField; pick a free temp value <= 32767
            _tmp = 32767
            while Isoform_Modification.objects.filter(
                isoform_id=_link.isoform_id, mod_position=_link.mod_position, mod_order=_tmp).exists():
                _tmp -= 1

            with transaction.atomic(using='dpeptide'):
                _link_order = _link.mod_order
                _other_order = _other.mod_order
                _link.mod_order = _tmp
                _link.save()
                _other.mod_order = _link_order
                _other.save()
                _link.mod_order = _other_order
                _link.save()
    return redirect('isoform_detail', pk=_isoform_pk)


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


#=================================================================================================
# Monomer
#=================================================================================================
class Monomer_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Monomer
    template_name = 'dpeptide/dict/dictionary_list.html'
    filterset_class = Monomer_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Monomer'
    app_name = 'dpeptide'
    ordering = ['symbol']

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Monomers'
        context['model_name'] = self.model_name
        context['list_url'] = 'monomer_list'
        context['create_url'] = 'monomer_create'
        context['update_url'] = 'monomer_update'
        return context


@login_required
def Monomer_CreateView(req):
    form = Monomer_Form(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user, str(instance.pk), 'Create Monomer', 'Completed')
                return redirect('monomer_list')
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('monomer_list')

    context = {'form': form, 'modal_title': 'Create Monomer', 'form_action': 'monomer_create'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


@login_required
def Monomer_UpdateView(req, pk):
    _object = get_object_or_404(Monomer, monomer_id=pk)
    if req.method == 'POST':
        form = Monomer_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dpeptide'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user, str(instance.pk), 'Update Monomer', 'Completed')
                return redirect('monomer_list')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('monomer_list')

    form = Monomer_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Monomer', 'modal_subtitle': pk, 'form_action': 'monomer_update'}
    return render(req, 'dpeptide/dpeptide_form_modal.html', context)


#=================================================================================================
# Sequence structure preview (BILN / HELM -> 2D)
#=================================================================================================
@login_required
def Sequence_StructureView(req):
    """Render a 2D structure SVG from a BILN or HELM string (via ?biln= / ?helm=)."""
    biln = req.GET.get('biln') or ''
    helm = req.GET.get('helm') or ''

    try:
        if biln:
            svg = render_biln_svg(biln)
        elif helm:
            svg = render_helm_svg(helm)
        else:
            svg = _error_svg('No BILN or HELM sequence provided')
        return HttpResponse(svg, content_type='image/svg+xml')
    except Exception as err:
        return HttpResponse(_error_svg(str(err)), content_type='image/svg+xml')
