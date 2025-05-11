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

# Create your views here.
from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from apputil.utils.filters_base import FilteredListView
from apputil.utils.views_base import permission_not_granted, HtmxupdateView, SimplecreateView, SimpleupdateView,  SimpledeleteView, CreateFileView

from adjcoadd.constants import *

#DScreen
from dscreen.models import  Screen_Run, Assay, AssayData_MIC, AssayData_CC50, AssayData_HC50
from dscreen.forms import (ScreenRun_Filter, 
                           #ScreenRun_DetailForm, ScreenRun_UpdateForm, ScreenRun_CreateForm,
                        )

#=================================================================================================
# Screen_Run
#=================================================================================================
class ScreenRun_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model = Screen_Run  
    template_name = 'dscreen/screenrun/screenrun_list.html'
    filterset_class = ScreenRun_Filter
    model_fields = model.HEADER_FIELDS
    model_name = 'Screen_Run'
    app_name = 'dscreen'
    ordering=['-acreated_at']

# -----------------------------------------------------------------
@login_required
def ScreenRun_DetailView(req, pk):
    context={}
    object_=get_object_or_404(Screen_Run, pk=pk)

#    smol_initial = Chem.MolToMolBlock(object_.smol) if object_.smol else None
#    form=Drug_form(instance=object_, initial={"smol":smol_initial},)
    form=ScreenRun_DetailForm(instance=object_,)
    context["object"]=object_
    context["form"]=form
    context["Links"]=LinkList

    #context['mol_img_url'] = settings.MOL_IMG_URL
    # try:
    #     context["object_mol"]=Chem.MolToMolBlock(object_.smol)
    #     m="\\n".join(context["object_mol"].split("\n"))
    #     context["object_mol"]=m
    # except Exception as err:
    #     context["object_mol"]=""
    return render(req, "ddrug/drug/drug_detail.html", context)


# -----------------------------------------------------------------
class ScreenRun_DeleteView(SimpledeleteView):
    model = Screen_Run
    transaction_use = 'dscreen'
