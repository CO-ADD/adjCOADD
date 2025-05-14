import os
import json
from rdkit import Chem
from django_filters.views import FilterView

from django.views.generic import ListView
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

# Create your views here.
from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from applib.django.views import Filtered_ListView, Base_CreateView, Base_UpdateView,  Base_DeleteView
#from applib.django.views import permission_not_granted, Htmx_UpdateView

from adjcoadd.constants import *

#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50
from dscreen.forms import (ScreenRun_Filter, ScreenRun_DetailForm,
                           #ScreenRun_UpdateForm, ScreenRun_CreateForm,
                        )

#=================================================================================================
# Screen_Run
#=================================================================================================
class ScreenRun_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model = Screen_Run  
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

# -----------------------------------------------------------------
@login_required
def ScreenRun_DetailView(req, pk):
    context={}
    object_=get_object_or_404(Screen_Run, pk=pk)

    form=ScreenRun_DetailForm(instance=object_,)
    context["object"]=object_
    context["form"]=form
    context["Links"]=LinkList

    return render(req, "dscreen/screenrun/screenrun_detail.html", context)

# -----------------------------------------------------------------
class ScreenRun_CreateView(Base_CreateView):
    form_class=ScreenRun_DetailForm
    template_name='dscreen/screenrun/screenrun_create.html'

# -----------------------------------------------------------------
class ScreenRun_UpdateView(Base_UpdateView):
    form_class=Screen_Run
    template_name='dscreen/screenrun/screenrun_update.html'
    model=Screen_Run

# -----------------------------------------------------------------
class ScreenRun_DeleteView(Base_DeleteView):
    model = Screen_Run
    transaction_use = 'dscreen'
