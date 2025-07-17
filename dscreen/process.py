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

# from adjcoadd.constants import *

from dscreen.models import Screen_Run

from dscreen.utils.screenrun_process import Upload_ReadOuts_Process, Upload_Motherplates_Process
from dscreen.utils.summary import update_screenrun_summary
from dsample.models import Project
from dplate.models import MasterPlate, TestPlate

from applib.process.process_forms import Process_View
from apputil.utils.form_wizard_tools import ImportHandler_View

from applib.process.process_forms import SelectSingleFile_StepForm,Finalize_StepForm,Upload_StepForm
#from apputil.utils.form_wizard_tools import SelectSingleFile_StepForm, Upload_StepForm, Finalize_StepForm 


class Readout_StepForm(SelectSingleFile_StepForm):
# --------------------------------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['multi_files'].label = 'Xlsx file from Tecan/BioTek readers'

class PlatePrep_StepForm(SelectSingleFile_StepForm):
# --------------------------------------------------------------------------------------------------
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['multi_files'].label = 'PlatePrep Xlsx Workbook'
        self.fields['multi_files'].help_text = mark_safe("Xlsx Workbook containing: <li> [TestPlateList] <li> [MotherPlates]")
        
# --------------------------------------------------------------------------------------------------
class Add_Readout_ProcessView(Process_View):
    process_name = 'Upload_Readout'
    model = Screen_Run

    #name_step1="Upload"
    form_list = [
        ('select_file', Readout_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]
    template_name = 'dscreen/screenrun_process/load_readouts.html'

    # customize util functions to validate files:
    # vitek -- upload_VitekPDF_Process
    def file_process_handler(self, request, *args, **kwargs):
        
        print(" [Add_Readouts.file_process_handler]")
        form_data=kwargs.get('form_data', None)
        
        _upload = False
        _overwrite=False
        
        valLog=Upload_ReadOuts_Process(request, self.file_dir, self.file_list, RunID=self.pk, upload=self.upload, appuser=request.user) 
        if self.upload:
            obj = Screen_Run.objects.select_for_update().get(run_id=self.pk)
            update_screenrun_summary(obj)
            obj.save(**kwargs)

        return(valLog)


# --------------------------------------------------------------------------------------------------
class Add_Motherplate_ProcessView(Process_View):
    process_name = 'Upload_Motherplate'
    model = Screen_Run

    form_list = [
        ('select_file', PlatePrep_StepForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dscreen/screenrun_process/load_motherplates.html'

    def file_process_handler(self, request, *args, **kwargs):    
        print(" [Add_Motherplate.file_process_handler]")
        form_data=kwargs.get('form_data', None)
        
        _upload = False
        _overwrite=False
        
        valLog=Upload_Motherplates_Process(request, self.file_dir, self.file_list, RunID=self.pk, upload=self.upload, appuser=request.user) 
 
        return(valLog)

# --------------------------------------------------------------------------------------------------
class Add_Testplate_ProcessView(Process_View):
    process_name = 'Upload_Testplates'
    model = Screen_Run
