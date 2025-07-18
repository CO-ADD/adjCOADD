import os
import json
import datetime

#from rdkit import Chem
#from django_filters.views import FilterView

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
from django.utils.safestring import mark_safe

from apputil.models import ApplicationLog
from apputil.forms import Document_Form
from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_RemoveView, Filtered_ListView
# Create your views here.

from dplate.models import TestPlate
from dplate.forms import TestPlate_Filter

#=================================================================================================
# TestPlates  
#=================================================================================================

class TestPlate_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model=TestPlate 
    template_name = 'dplate/testplate/testplate_list.html' 
    filterset_class=TestPlate_Filter
    model_fields=model.LIST_VIEW_FIELDS
    model_name = 'TestPlate'
    app_name = 'dplate'

#-------------------------------------------------------------------------------------------------
def TestPlate_MapView(req,pk):
    print(f" [Testplate] MapView {pk}  {req}")
    
    context = {}
    if req.method=='GET':
        _object=get_object_or_404(TestPlate, plate_id=pk)
        _now = datetime.datetime.now()
        _xls_name = f'Testplate_{pk}_{_now:%Y%m%d}.xlsx'
        
        return render(req,'dplate/testplate/testplate_map.html',context)
