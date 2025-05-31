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

from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_DeleteView, Filtered_ListView
from apputil.utils.form_wizard_tools import ImportHandler_View, SelectMultipleFiles_StepForm, Upload_StepForm, Finalize_StepForm

# from adjcoadd.constants import *

from dscreen.models import Screen_Run
from dsummary.models import Summary_ScreenRun
from dscreen.forms import ScreenRun_Filter, ScreenRun_CreateForm, ScreenRun_UpdateForm
from dsample.models import Project
from dplate.models import MasterPlate, TestPlate

#=================================================================================================
# ScreenRun
#=================================================================================================
class ScreenRun_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Summary_ScreenRun  
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

