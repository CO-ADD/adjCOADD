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

from dcollab.models import Organisation, Collab_Group, Collab_User
from dcollab.forms import Organisation_Filter, Project_CreateForm, Project_UpdateForm
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

