import os
import json
import datetime

#from rdkit import Chem
#from django_filters.views import FilterView

# from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
# from django.contrib.auth.mixins import LoginRequiredMixin
# from django.contrib import messages
# from django.core.exceptions import ValidationError
# from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
# from django.db import transaction, IntegrityError
# from django.db.models import Count
# from django.http import JsonResponse
# from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
# from django.urls import reverse_lazy
# from django.utils.functional import SimpleLazyObject
from django.utils.safestring import mark_safe
from django import forms

# from apputil.models import ApplicationLog
# from apputil.forms import Document_Form
# from applib.django.base.views import Base_CreateView, Base_UpdateView, Base_RemoveView, Filtered_ListView

from adjcoadd.constants import *
from dscreen.models import Screen_Run
from dscreen.utils.screenrun_process import (Summary_ScreenRun_Process, 
                                            Upload_ReadOuts_Process, Upload_Motherplates_Process,Upload_TestplateList_Process)
# from dsample.models import Project
# from dplate.models import MasterPlate, TestPlate

from applib.process.process_forms import Process_View
from apputil.utils.form_wizard_tools import ImportHandler_View

from applib.process.process_forms import SelectSingleFile_StepForm, Finalize_StepForm, Upload_StepForm
#from apputil.utils.form_wizard_tools import SelectSingleFile_StepForm, Upload_StepForm, Finalize_StepForm 


# --------------------------------------------------------------------------------------------------
class PlatePrep_SelectForm(SelectSingleFile_StepForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['multi_files'].label = 'PlatePrep Workbook'

        _help_text = 'XLSX Workbook containing the following Sheets: '
        _help_text += '<li> [MotherPlates]'
        _help_text += '<li> [TestPlateList, Assays, (HCPrep, PSPrep)]'
        # _help_text += '<p> Use <a href "{% static '
        # _help_text += f"'{DOC_TEMPLATES['hc_plateprep']}'"
        # _help_text += ' %}"> PlatePrep Template</a> '
        #_help_text += f'<p> Use <a href="static/{DOC_TEMPLATES["hc_plateprep"]}">PlatePrep Template</a>'
        #_help_text += ' {% static '
        #_help_text += f"'{DOC_TEMPLATES['hc_plateprep']}'"
        
        self.fields['multi_files'].help_text = mark_safe(_help_text)

        #<a href="{% static 'django-pdf/generator/static/pdfs/nowy.pdf' %}">{{ file }}</a>

# --------------------------------------------------------------------------------------------------
class Readout_SelectForm(SelectSingleFile_StepForm):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['multi_files'].label = 'Xlsx file from Tecan/BioTek readers'


# --------------------------------------------------------------------------------------------------
class Load_Readouts_ProcessView(Process_View):
    process_name = 'Upload_Readouts'
    model = Screen_Run

    #name_step1="Upload"
    form_list = [
        ('select_file', Readout_SelectForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dscreen/screenrun_process/load_readouts.html'

    select_html  = 'Please select a Excel [xlsx] file from Tecan/BioTek readers'
    select_html += '\n Make sure file contains correct  <b>TestPlate IDs</b>'

    upload_html  = 'Please check the TestPlate IDs [<i>Item</i>] for any "New Testplate" [<i>Action</i>]'
    upload_html += '\n Make sure the IDs are unique and reflect the IDs in <b>TestPLateList</b>'
    upload_html += '\n In case, correct the IDs in the <b>Readout</b> file and repeat the upload'

    message_html =[
       ('select_file',select_html),
       ('upload',upload_html),
       ('finalize','') 
    ]

    # Customize Function to validate and upload files:
    def file_process_handler(self, request, *args, **kwargs):

        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'overwrite' in form_data:
            self.overwrite = form_data['overwrite']

        valLog=Upload_ReadOuts_Process(request, self.file_dir, self.file_list, RunID=self.pk, 
                                       upload=self.upload, overwrite=self.overwrite,appuser=request.user)

        return(valLog)

    # Customize Function to update after upload:
    def file_process_finalizer(self, request, pk):
        Summary_ScreenRun_Process(request, pk)

# --------------------------------------------------------------------------------------------------
class Load_Motherplates_ProcessView(Process_View):
    process_name = 'Upload_Motherplates'
    model = Screen_Run

    form_list = [
        ('select_file', PlatePrep_SelectForm),
        ('upload', Upload_StepForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dscreen/screenrun_process/load_motherplates.html'

    def file_process_handler(self, request, *args, **kwargs):    
        print(" [Add_Motherplate.file_process_handler]")

        self.upload = False
        self.overwrite = False
        
        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'overwrite' in form_data:
            self.overwrite = form_data['overwrite']
       
        valLog=Upload_Motherplates_Process(request, self.file_dir, self.file_list, RunID=self.pk, 
                                           upload=self.upload, overwrite=self.overwrite, appuser=request.user) 
 
        return(valLog)

    # Customize Function to update after upload:
    def file_process_finalizer(self, request, pk):
        Summary_ScreenRun_Process(request, pk)

# --------------------------------------------------------------------------------------------------
class TestPlate_UploadForm(forms.Form):
    apply_mp = forms.BooleanField(initial=False, required=False, help_text="Fill Testplates with Compounds, Layout and Assays")
#    only_dr = forms.BooleanField(initial=False, required=False, help_text="Recalculate only Doseresponse (after change of PlateQuality)")
    upload = forms.BooleanField(initial=False, required=False, help_text="Upload Testplates and Doseresponse to Database")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['apply_mp'].label = "Apply Motherplates/Layout"
        # self.fields['only_dr'].label = "Only Doseresponse"
        self.fields['upload'].label = "Upload Data"


# --------------------------------------------------------------------------------------------------
class Load_TestplateList_ProcessView(Process_View):
    process_name = 'Upload_TestplatesList'
    model = Screen_Run
    progress = [0,0]

    form_list = [
        ('select_file', PlatePrep_SelectForm),
        ('upload', TestPlate_UploadForm),
        ('finalize', Finalize_StepForm),
    ]

    template_name = 'dscreen/screenrun_process/load_testplatelist.html'

    # select_html  = 'Please select a Excel [xlsx] file from Tecan/BioTek readers'
    # select_html += '\n Make sure file contains correct  <b>TestPlate IDs</b>'

    # upload_html  = 'Please check the TestPlate IDs [<i>Item</i>] for any "New Testplate" [<i>Action</i>]'
    # upload_html += '\n Make sure the IDs are unique and reflect the IDs in <b>TestPLateList</b>'
    # upload_html += '\n In case, correct the IDs in the <b>Readout</b> file and repeat the upload'

    # message_html =[
    #    ('select_file',select_html),
    #    ('upload',upload_html),
    #    ('finalize','') 
    # ]
    def file_process_handler(self, request, *args, **kwargs):    
        #print(" [Add_TestplateList_ProcessView.file_process_handler]")

        self.upload = False
        self.apply_mp = False
        # self.only_dr = False
                
        # Set Form Data        
        form_data=kwargs.get('form_data', None)
        if 'upload' in form_data:
            self.upload = form_data['upload']
        if 'apply_mp' in form_data:
            self.apply_mp = form_data['apply_mp']
        # if 'only_dr' in form_data:
        #     self.only_dr = form_data['only_dr']

        valLog=Upload_TestplateList_Process(request, self.file_dir, self.file_list, RunID=self.pk, 
                                           upload=self.upload, apply_mp=self.apply_mp, 
                                           appuser=request.user) 
 
        return(valLog)

    def file_process_finalizer(self, request, pk):
        Summary_ScreenRun_Process(request, pk)
