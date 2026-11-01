import datetime
import io

import numpy as np
import matplotlib.pyplot as plt

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import FieldDoesNotExist
from django.db import transaction, IntegrityError
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, render, redirect

from apputil.models import ApplicationLog, Dictionary
from applib.django.base.views import Filtered_ListView

from dplate.models import TestPlate, MasterPlate, Labware
from dplate.forms import (TestPlate_Filter, MasterPlate_Filter, Labware_Filter,
                          TestPlate_UpdateForm, MasterPlate_UpdateForm,
                          Labware_Form, Labware_UpdateForm)


#=================================================================================================
# TestPlate
#=================================================================================================

TestPlate_DETAIL_FIELDS = ['plate_id', 'plate_type', 'labware_id', 'n_rows', 'n_cols', 'n_wells',
                           'run_id', 'assay_id', 'result_type', 'plate_quality', 'test_date',
                           'control_layout', 'zfactor', 'plate_qc', 'process_status',
                           'n_reads', 'n_samples', 'n_layout', 'n_inhibitions',
                           'n_doseresponses', 'n_synergies']


class TestPlate_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = TestPlate
    template_name = 'dplate/plate_list.html'
    filterset_class = TestPlate_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'TestPlate'
    app_name = 'dplate'
    ordering = ['plate_id']
    paginate_by = 50

    def get_paginate_by(self, queryset):
        if self.request.GET.get('view') == 'card':
            return 15
        return super().get_paginate_by(queryset)

    def get_queryset(self):
        # Filter-first: the testplate table is very large.
        if not self.get_filter_record():
            self.filterset = self.filterset_class(self.request.GET,
                                                  queryset=self.model.objects.none(),
                                                  filterset_dict={}, deep=False)
            self.filter_count = 0
            return self.model.objects.none()
        return super().get_queryset()

    def get_context_data(self, **kwargs):
        context = super(Filtered_ListView, self).get_context_data(**kwargs)
        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        context['filter'] = self.filterset
        context['paginate_by'] = self.get_paginate_by(self.object_list)
        context['fields'] = self.model.get_fields(fields=self.model_fields)
        context['filterset'] = filter_record
        context['querycount'] = self.filter_count
        context['Count'] = None
        context['has_filters'] = bool(filter_record_dict)
        context['model_name'] = self.model_name
        context['view'] = self.request.GET.get('view', 'table')
        context['list_url'] = 'testplate_list'
        context['detail_url'] = 'testplate_detail'
        return context


# -----------------------------------------------------------------
def TestPlate_MapView(req, pk):
    _object = get_object_or_404(TestPlate, plate_id=pk)
    context = {'object': _object}
    return render(req, 'dplate/testplate/testplate_map.html', context)


# -----------------------------------------------------------------
def _plate_map_data(instance):
    """Build a JSON-safe payload describing a plate and its wells for plate_map.js."""
    is_test = isinstance(instance, TestPlate)
    instance.get_wells()
    row_labels = instance.ROW_LABELS[:instance.n_rows]
    wells = []

    for well_id in sorted(instance.wells.keys()):
        well = instance.wells.get(well_id)
        if well is None:
            continue
        _r, _c = instance.well_rowcol(well_id)

        entry = {
            'id': well_id,
            'row': row_labels[_r - 1],
            'col': _c,
            'cmpbatch': [str(x) for x in (well.cmpbatch_lst or []) if x not in ('', None)],
            'conc': [float(x) for x in (well.conc_lst or []) if x not in ('', None)],
            'conc_unit': [str(x) for x in (well.conc_unit_lst or []) if x != ''],
        }

        if is_test:
            entry['readout'] = [float(x) for x in (well.readouts or []) if x not in ('', None)]
            entry['readout_types'] = [str(x) for x in (well.readout_types or [])]
            entry['inhibition'] = float(well.inhibition) if well.inhibition not in (None, -1) else None
            entry['act_type'] = well.act_type or ''
            entry['is_control'] = bool(well.is_control)
            entry['is_poscontrol'] = bool(well.is_poscontrol)
            entry['is_negcontrol'] = bool(well.is_negcontrol)
            entry['is_sample'] = bool(well.is_sample)
            entry['is_skip'] = bool(well.is_skip)
        else:
            entry['barcode'] = well.barcode or ''

        wells.append(entry)

    metrics = {}
    if is_test:
        max_readouts = 0
        for w in wells:
            if len(w['readout']) > max_readouts:
                max_readouts = len(w['readout'])
        for i in range(max_readouts):
            type_label = ''
            for w in wells:
                rts = w.get('readout_types') or []
                if i < len(rts) and rts[i]:
                    type_label = rts[i]
                    break
            metrics[f'readout{i}'] = {
                'label': f'Readout {i + 1}',
                'type_label': type_label,
                'mode': 'seq', 'key': 'readout', 'index': i,
            }
        metrics['inhib'] = {'label': 'Inhibition (%)', 'mode': 'diverging', 'min': 0, 'max': 100, 'key': 'inhibition'}

    # Highest concentration index (1-based) with a non-zero value across wells
    max_conc = 0
    for w in wells:
        for i, c in enumerate(w['conc']):
            if c and i + 1 > max_conc:
                max_conc = i + 1
    for i in range(max_conc):
        metrics[f'conc{i}'] = {'label': f'Concentration {i + 1}', 'mode': 'seq', 'key': 'conc', 'index': i}

    return {
        'plate_id': instance.plate_id,
        'n_rows': instance.n_rows,
        'n_cols': instance.n_cols,
        'n_wells': instance.n_wells,
        'row_labels': row_labels,
        'show_type': is_test,
        'metrics': metrics,
        'wells': wells,
    }


# -----------------------------------------------------------------
def _stat(arr, idx):
    """Safely return an ArrayField element (e.g. median/MAD) for display."""
    if arr and len(arr) > idx:
        val = arr[idx]
        if val is not None:
            return val
    return '—'


def _quality_class(value):
    """Traffic-light CSS class for a plate quality dictionary value."""
    v = (value or '').lower()
    if 'reject' in v:
        return 'q-rejected'
    if 'retest' in v:
        return 'q-retest'
    if 'invalid' in v:
        return 'q-invalid'
    if 'valid' in v:
        return 'q-valid'
    return ''


# -----------------------------------------------------------------
def _build_display_groups(instance, field_names, chunk=6):
    groups = []
    rows = []
    for name in field_names:
        try:
            field = instance._meta.get_field(name)
        except FieldDoesNotExist:
            continue
        value = getattr(instance, name, None)
        if value is None or value == '':
            value = ''
        else:
            value = str(value)
        rows.append({'label': field.verbose_name, 'value': value})
    for i in range(0, len(rows), chunk):
        groups.append(rows[i:i + chunk])
    return groups


@login_required
def TestPlate_DetailView(req, pk):
    _object = get_object_or_404(TestPlate, plate_id=pk)

    # Previous/Next TestPlate within the same Screen_Run (same ordering as ScreenRun_DetailView)
    prev_plate = None
    next_plate = None
    if _object.run_id:
        _siblings = list(
            TestPlate.objects.filter(run_id=_object.run_id, astatus__gt=-1)
            .order_by('plate_id')
            .values_list('plate_id', flat=True)
        )
        if _object.plate_id in _siblings:
            _idx = _siblings.index(_object.plate_id)
            if _idx > 0:
                prev_plate = _siblings[_idx - 1]
            if _idx < len(_siblings) - 1:
                next_plate = _siblings[_idx + 1]

    quality_choices = [
        (d.dict_value, d.strtml(), _quality_class(d.dict_value))
        for d in Dictionary.get_filterobj(TestPlate.DICTIONARY_FIELDS['plate_quality'])
    ]
    _MED = TestPlate.STATS_MEDIAN
    _MAD = TestPlate.STATS_MAD
    plate_info = {
        'assay_id': str(_object.assay_id) if _object.assay_id else '—',
        'zfactor': _object.zfactor if _object.zfactor not in (None, -1) else '—',
        'poscontrol': f"{_stat(_object.poscontrol_stats, _MED)} +/- {_stat(_object.poscontrol_stats, _MAD)}",
        'negcontrol': f"{_stat(_object.negcontrol_stats, _MED)} +/- {_stat(_object.negcontrol_stats, _MAD)}",
        'sample': f"{_stat(_object.sample_stats, _MED)} +/- {_stat(_object.sample_stats, _MAD)}",
        'edge': f"Edge: {_stat(_object.edge_stats, 0)} - NonEdge: {_stat(_object.edge_stats, 1)}",
        'plate_quality': _object.plate_quality_id or '',
        'quality_choices': quality_choices,
    }
    context = {
        'object': _object,
        'display_groups': _build_display_groups(_object, TestPlate_DETAIL_FIELDS),
        'list_url': 'testplate_list',
        'update_url': 'testplate_update',
        'plate_map_data': _plate_map_data(_object),
        'plate_info': plate_info,
        'prev_plate': prev_plate,
        'next_plate': next_plate,
    }
    return render(req, 'dplate/plate_detail.html', context)


@login_required
def TestPlate_UpdateView(req, pk):
    _object = get_object_or_404(TestPlate, plate_id=pk)
    if req.method == 'POST':
        form = TestPlate_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dplate'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user,
                                       str(instance.pk), 'Update TestPlate', 'Completed')
                return redirect('testplate_detail', pk=instance.pk)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('testplate_detail', pk=pk)

    form = TestPlate_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update TestPlate', 'modal_subtitle': pk}
    return render(req, 'dplate/plate_form_modal.html', context)


@login_required
def TestPlate_QualityUpdateView(req, pk):
    _object = get_object_or_404(TestPlate, plate_id=pk)
    if req.method == 'POST':
        quality_value = (req.POST.get('plate_quality') or '').strip()
        if quality_value:
            _dict = Dictionary.get(TestPlate.DICTIONARY_FIELDS['plate_quality'], quality_value)
            if _dict:
                try:
                    with transaction.atomic(using='dplate'):
                        _object.plate_quality = _dict
                        _object.save(user=req.user)
                        ApplicationLog.add('Update', str(_object.pk), 'Info', req.user,
                                           str(_object.pk), f'Set TestPlate quality {quality_value}', 'Completed')
                    messages.success(req, f'Plate quality set to {quality_value}')
                except Exception as err:
                    messages.warning(req, f'Update failed due to {err} error')
            else:
                messages.warning(req, f'Unknown plate quality: {quality_value}')
        else:
            messages.warning(req, 'No plate quality provided')
    return redirect('testplate_detail', pk=pk)


def _empty_heatmap_png():
    """Return a transparent 1px PNG for plates with no readout data."""
    buf = io.BytesIO()
    fig, ax = plt.subplots(figsize=(1, 1), dpi=50)
    ax.axis('off')
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    fig.savefig(buf, format='png', transparent=True)
    plt.close(fig)
    buf.seek(0)
    return HttpResponse(buf.read(), content_type='image/png')


@login_required
def TestPlate_HeatmapView(req, pk):
    """Static mini heatmap of readout[0] for TestPlate card views (no axes labels)."""
    _object = get_object_or_404(TestPlate, plate_id=pk)
    _object.get_wells()

    rows = _object.n_rows or 0
    cols = _object.n_cols or 0
    if rows <= 0 or cols <= 0:
        return _empty_heatmap_png()

    mat = np.full((rows, cols), np.nan)
    has_value = False
    for well_id, well in _object.wells.items():
        if well is None:
            continue
        try:
            r, c = _object.well_rowcol(well_id)
            if well.readouts:
                mat[r - 1, c - 1] = float(well.readouts[0])
                has_value = True
        except Exception:
            continue

    if not has_value:
        return _empty_heatmap_png()

    # readout type label: plate-level first, else first well's readout_types[0]
    _rt = _object.readout_type or ''
    if not _rt:
        for well in _object.wells.values():
            if well is not None and well.readout_types and well.readout_types[0]:
                _rt = well.readout_types[0]
                break

    fig_w = max(2.5, cols * 0.18)
    fig_h = max(2.0, rows * 0.18) + 0.3
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=100)
    cmap = plt.cm.Reds.copy()
    cmap.set_bad('#f0f0f0')
    vmin = float(np.nanmin(mat))
    vmax = float(np.nanmax(mat))
    if vmin == vmax:
        vmax = vmin + 1
    ax.imshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, aspect='equal', origin='upper')
    ax.set_title(_rt or 'Readout', fontsize=9, pad=2)
    ax.axis('off')
    fig.subplots_adjust(left=0.02, right=0.98, top=0.84, bottom=0.02)

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=100)
    plt.close(fig)
    buf.seek(0)
    return HttpResponse(buf.read(), content_type='image/png')


#=================================================================================================
# MasterPlate
#=================================================================================================

MasterPlate_DETAIL_FIELDS = ['plate_id', 'plate_type', 'labware_id', 'n_rows', 'n_cols', 'n_wells',
                             'run_id', 'prep_date', 'plate_quality', 'plating', 'cpoz_id', 'well_type',
                             'dilution_layout', 'master_solvent', 'master_conc', 'master_conc_unit',
                             'master_volume', 'master_volume_unit', 'n_samples', 'n_layout', 'process_status']


class MasterPlate_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = MasterPlate
    template_name = 'dplate/plate_list.html'
    filterset_class = MasterPlate_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'MasterPlate'
    app_name = 'dplate'
    ordering = ['plate_id']
    paginate_by = 50

    def get_paginate_by(self, queryset):
        if self.request.GET.get('view') == 'card':
            return 15
        return super().get_paginate_by(queryset)

    def get_queryset(self):
        # Filter-first: the masterplate table is very large.
        if not self.get_filter_record():
            self.filterset = self.filterset_class(self.request.GET,
                                                  queryset=self.model.objects.none(),
                                                  filterset_dict={}, deep=False)
            self.filter_count = 0
            return self.model.objects.none()
        return super().get_queryset()

    def get_context_data(self, **kwargs):
        context = super(Filtered_ListView, self).get_context_data(**kwargs)
        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: " + str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        context['filter'] = self.filterset
        context['paginate_by'] = self.get_paginate_by(self.object_list)
        context['fields'] = self.model.get_fields(fields=self.model_fields)
        context['filterset'] = filter_record
        context['querycount'] = self.filter_count
        context['Count'] = None
        context['has_filters'] = bool(filter_record_dict)
        context['model_name'] = self.model_name
        context['view'] = self.request.GET.get('view', 'table')
        context['list_url'] = 'masterplate_list'
        context['detail_url'] = 'masterplate_detail'
        return context


@login_required
def MasterPlate_DetailView(req, pk):
    _object = get_object_or_404(MasterPlate, plate_id=pk)
    context = {
        'object': _object,
        'display_groups': _build_display_groups(_object, MasterPlate_DETAIL_FIELDS),
        'list_url': 'masterplate_list',
        'update_url': 'masterplate_update',
        'plate_map_data': _plate_map_data(_object),
    }
    return render(req, 'dplate/plate_detail.html', context)


@login_required
def MasterPlate_UpdateView(req, pk):
    _object = get_object_or_404(MasterPlate, plate_id=pk)
    if req.method == 'POST':
        form = MasterPlate_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dplate'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user,
                                       str(instance.pk), 'Update MasterPlate', 'Completed')
                return redirect('masterplate_detail', pk=instance.pk)
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('masterplate_detail', pk=pk)

    form = MasterPlate_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update MasterPlate', 'modal_subtitle': pk}
    return render(req, 'dplate/plate_form_modal.html', context)


#=================================================================================================
# Labware
#=================================================================================================

class Labware_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Labware
    template_name = 'dplate/labware/labware_list.html'
    filterset_class = Labware_Filter
    model_fields = model.LIST_VIEW_FIELDS
    model_name = 'Labware'
    app_name = 'dplate'
    ordering = ['labware_id']
    paginate_by = 50

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['model_name'] = self.model_name
        context['list_url'] = 'labware_list'
        context['create_url'] = 'labware_create'
        context['update_url'] = 'labware_update'
        return context


@login_required
def Labware_CreateView(req):
    form = Labware_Form(req.POST or None)
    if req.method == 'POST':
        if form.is_valid():
            try:
                with transaction.atomic(using='dplate'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Create', str(instance.pk), 'Info', req.user,
                                       str(instance.pk), 'Create Labware', 'Completed')
                return redirect('labware_list')
            except IntegrityError as err:
                messages.error(req, f'IntegrityError {err} happens, record may be existed!')
        else:
            messages.warning(req, f'Create failed due to {form.errors} error')
        return redirect('labware_list')

    context = {'form': form, 'modal_title': 'Create Labware'}
    return render(req, 'dplate/plate_form_modal.html', context)


@login_required
def Labware_UpdateView(req, pk):
    _object = get_object_or_404(Labware, labware_id=pk)
    if req.method == 'POST':
        form = Labware_UpdateForm(req.POST, instance=_object)
        if form.is_valid():
            try:
                with transaction.atomic(using='dplate'):
                    instance = form.save(commit=False)
                    instance.save(user=req.user)
                    ApplicationLog.add('Update', str(instance.pk), 'Info', req.user,
                                       str(instance.pk), 'Update Labware', 'Completed')
                return redirect('labware_list')
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
        else:
            messages.warning(req, f'Update failed due to {form.errors} error')
        return redirect('labware_list')

    form = Labware_UpdateForm(instance=_object)
    context = {'object': _object, 'form': form, 'modal_title': 'Update Labware', 'modal_subtitle': pk}
    return render(req, 'dplate/plate_form_modal.html', context)

