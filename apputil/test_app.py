import os
from formtools.wizard.views import SessionWizardView
from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse, QueryDict
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView
from django.views.generic.detail import DetailView
from django.db import transaction, IntegrityError

# from adjcoadd.constants import *
# from dorganism.models import Organism, Taxonomy
# from ddrug.models import Drug, VITEK_Card, VITEK_AST, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
# from dgene.models import Genome_Sequence,Gene, WGS_CheckM, WGS_FastQC, ID_Pub, ID_Sequence
# from dscreen.models import Screen_Run
# from dsample.models import Project
# from dcell.models import Cell

# from apputil.forms import Login_Form, AppUser_Form, AppUser_Filter, AppLog_Filter, Dictionary_Filter, Dictionary_Form, Document_Form 
# from apputil.models import ApplicationUser, Dictionary, ApplicationLog, Document
# from applib.django.base.views import (SuperUserRequiredMixin, permission_not_granted, 
#                                       Base_CreateView, Base_UpdateView, Base_RemoveView,
#                                       Filtered_ListView, 
#                                       Htmx_UpdateView, File_CreateView, Base_DataExportView)

# from apputil.utils.files_upload import Importhandler, OverwriteStorage, file_location
# from apputil.utils.data_style import convert_heatmap
# from apputil.utils.form_wizard_tools import SelectMultipleFiles_StepForm,SelectSingleFile_StepForm
# from apputil.utils.validation_log import Validation_Log

#-------------------------------------------------------------------------------------------------
def Test_View(req):
    context = {}
    if req.method=='GET':
        return render(req,'test/test_modal.html',context)

    
