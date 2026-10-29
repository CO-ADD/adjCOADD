import datetime

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import FieldDoesNotExist
from django.db import transaction
from django.shortcuts import get_object_or_404, render, redirect

from apputil.models import ApplicationLog
from applib.django.base.views import Filtered_ListView

from dplate.models import TestPlate, MasterPlate, Labware
from dplate.forms import (TestPlate_Filter, MasterPlate_Filter, Labware_Filter,
                          TestPlate_UpdateForm, MasterPlate_UpdateForm, Labware_UpdateForm)


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
    context = {
        'object': _object,
        'display_groups': _build_display_groups(_object, TestPlate_DETAIL_FIELDS),
        'list_url': 'testplate_list',
        'update_url': 'testplate_update',
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
        context['update_url'] = 'labware_update'
        return context


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

