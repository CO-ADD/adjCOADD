'''
View for uploading Vitek PDFs
'''
import pandas as pd

from django import forms
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.http import JsonResponse
from django.shortcuts import HttpResponse, render, redirect

from apputil.utils.files_upload import file_location, OverwriteStorage
from apputil.utils.form_wizard_tools import ImportHandler_WizardView, UploadFileForm, StepForm_1, FinalizeForm
from apputil.utils.validation_log import * 
from dorganism.models import Organism_Batch
from ddrug.models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub
from ddrug.utils.vitek import upload_VitekPDF_List

#==  VITEK Import View =============================================================

# customized Form
class VitekValidation_StepForm(StepForm_1):
    orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',
                                       widget=forms.Select(attrs={'class': 'form-control'}), 
                                       required=False, 
                                       help_text='(**Optional, if choose an Organism Batch ID, it will be applied in all uploaded files)',
                                       queryset=Organism_Batch.objects.filter(astatus__gte=0), )
    field_order = ['orgbatch_id', 'confirm']
# 

# Process View
class Import_VitekView(ImportHandler_WizardView):
    
    name_step1="Validation" # step label in template
    # define more steps name
    #... 
    # define each step's form
    form_list = [
        ('upload_file', UploadFileForm),
        ('step1', VitekValidation_StepForm),
        # add more step -> StepForm
        ('finalize', FinalizeForm),
    ]
    # define template
    template_name = 'ddrug/importhandler_vitek.html'
    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filelist=[]
        self.orgbatch_id=None
        self.valLog=None
        self.upload=False

        
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request
        # config session save and process cancel 
        SessionKey=self.request.session.session_key
        if current_step == 'upload_file':
            self.storage.extra_data['validation_result']="***have no result, please uploading files***"
            DirName = file_location(instance=request.user)  # define file store path during file process
            files = []
            if form.is_valid():
                if 'upload_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('upload_file-multi_files'))          
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=DirName)
                    filename = fs.save(f.name, f)
                    self.filelist.append(filename)

                # Parse PDF and Validation
                self.valLog=upload_VitekPDF_List(request, DirName, self.filelist, SessionKey=SessionKey, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user) 
                if self.valLog.nLogs['Error'] >0 :
                    dfLog = pd.DataFrame(self.valLog.get_aslist(logTypes= ['Error']))#convert result in a table
                    self.storage.extra_data['Confirm_to_Save'] = False
                else:
                    dfLog = pd.DataFrame(self.valLog.get_aslist())
                    self.storage.extra_data['Confirm_to_Save'] = True
                html=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False)
                self.storage.extra_data['validation_result']=html      
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['DirName'] = DirName          
            else:
                return render(request, 'ddrug/importhandler_vitek.html', context)

        elif current_step == 'step1': # recheck and save to DB
            # print("step validation again")
            self.upload=self.storage.extra_data['Confirm_to_Save']
            form = VitekValidation_StepForm(request.POST)
            self.organism_batch=request.POST.get("upload_file-orgbatch_id") #get organism_batch
            DirName=self.storage.extra_data['DirName'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name   
       
            self.valLog=upload_VitekPDF_List(request, DirName, self.filelist, SessionKey=SessionKey, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user)
            if self.valLog.nLogs['Error'] >0 :
                dfLog = pd.DataFrame(self.valLog.get_aslist(logTypes= ['Error']))#convert result in a table
            else:
                dfLog = pd.DataFrame(self.valLog.get_aslist())
            html=dfLog.to_html(classes=[ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"], index=False)
            self.storage.extra_data['validation_result']=html  
                            
        return self.get_form_step_data(form)

    def done(self, form_list, **kwargs):
        print("Finalize")
        # Redirect to the desired page after finishing
        filelist=self.storage.extra_data['filelist']
        for f in filelist:
            self.delete_file(f)
            print(filelist)
        return redirect(self.request.META['HTTP_REFERER'])  

    def get_context_data(self, form, **kwargs):
        context = super().get_context_data(form=form, **kwargs)
        # save information to context,
        # then display in templates  
        context['step1']=self.name_step1
        current_step = self.steps.current  
        if current_step == 'upload_file':
            context['validation_result']="***have no result, please finish uploading files firstly***"
        else:
            context['validation_result'] = self.storage.extra_data.get('validation_result', None)# 
        # if current_step == 'step1':
            context['Confirm_to_Save']=self.storage.extra_data.get('Confirm_to_Save', None)
        return context
    