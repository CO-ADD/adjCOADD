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
from apputil.utils.form_wizard_tools import ImportHandler_WizardView, SelectFile_StepForm, Upload_StepForm, StepForm_1, Finalize_StepForm
from apputil.utils.validation_log import * 
from dorganism.models import Organism_Batch
from ddrug.models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub
from ddrug.utils.vitek import upload_VitekPDF_Process

#==  VITEK Import Wizard =============================================================

# customized Form
class VitekValidation_StepForm(Upload_StepForm):
    orgbatch_id=forms.ModelChoiceField(label='Choose an Organism Batch',
                                       widget=forms.Select(attrs={'class': 'form-control'}), 
                                       required=False, 
                                       help_text='(**Optional, if choose an Organism Batch ID, it will be applied in all uploaded files)',
                                       queryset=Organism_Batch.objects.filter(astatus__gte=0), )
    field_order = ['orgbatch_id', 'confirm']
# 

# Process View
class Import_VitekView(ImportHandler_WizardView):
    
    name_step1="Upload" # step label in template
    # define more steps name
    #... 
    # define each step's form
    form_list = [
        ('select_file', SelectFile_StepForm),
        ('upload', VitekValidation_StepForm),
        # add more step -> StepForm
        ('finalize', Finalize_StepForm),
    ]
    # define template
    template_name = 'ddrug/importhandler_vitek.html'
    # Define a file storage for handling file uploads
    file_storage = FileSystemStorage(location='/tmp/')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filelist=[]
        self.dirname=None
        self.orgbatch_id=None
        self.valLog=None
        self.upload=False
        self.storage = None
        self.html_columns = ['Type','Note','Item','Filename','Help']
        # self.html_classes = [ "table", "table-bordered", "fixTableHead", "bg-light", "m-0"]
         
    def process_step(self, form):
        current_step = self.steps.current
        request = self.request
        # config session save and process cancel 
        SessionKey=self.request.session.session_key
        if current_step == 'select_file':
            self.storage.extra_data['validation_result']="-"
            self.dirname = file_location(instance=request.user)  # define file store path during file process
            files = []
            if form.is_valid():
                if 'select_file-multi_files' in request.FILES:
                    files.extend(request.FILES.getlist('select_file-multi_files'))          
                # Get clean FileList
                for f in files:
                    fs = OverwriteStorage(location=self.dirname)
                    filename = fs.save(f.name, f)
                    self.filelist.append(filename)

                # Parse PDF and Validation
                self.valLog=upload_VitekPDF_Process(request, self.dirname, self.filelist, SessionKey=SessionKey, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user) 
                if self.valLog.nLogs['Error'] >0 :
                    dfLog = self.valLog.get_ashtml(logTypes= ['Error'], columns=self.html_columns)#convert result in a table
                    self.storage.extra_data['confirm_to_upload'] = False
                else:
                    dfLog = self.valLog.get_ashtml(columns=self.html_columns)
                    self.storage.extra_data['confirm_to_upload'] = True

                self.storage.extra_data['validation_result'] = dfLog
                self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) checked for errors." 
                self.storage.extra_data['filelist'] = self.filelist
                self.storage.extra_data['dirname'] = self.dirname          
            else:
                self.storage.extra_data['validation_result']="No files selected"
                return render(request, 'ddrug/importhandler_vitek.html', context)

        elif current_step == 'upload': # recheck and save to DB
            # print("step validation again")
            form = VitekValidation_StepForm(request.POST)
            self.organism_batch=request.POST.get("upload_file-orgbatch_id") #get organism_batch
            self.upload=True
            self.dirname=self.storage.extra_data['dirname'] #get file path
            self.filelist=self.storage.extra_data['filelist'] #get files' name   
       
            self.valLog=upload_VitekPDF_Process(request, self.dirname, self.filelist, SessionKey=SessionKey, OrgBatchID=self.orgbatch_id, upload=self.upload, appuser=request.user)

            if self.valLog.nLogs['Error'] >0 :
                dfLog = self.valLog.get_asdf(logTypes= ['Error']) 
            else:
                dfLog = self.valLog.get_asdf()

            html=dfLog.to_html(columns=self.html_columns,classes=self.html_classes,index=False).replace("\\n","<br>")
            self.storage.extra_data['validation_result'] = html  
            self.storage.extra_data['validation_message']= f" {len(self.filelist)} file(s) Uploaded." 
                            
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
        context['validation_message'] = self.storage.extra_data.get('validation_message', None)
        if current_step == 'upload_file':
            context['validation_result']="Select VITEK PDF files"
        else:
            context['validation_result'] = self.storage.extra_data.get('validation_result', None) 
        # if current_step == 'step1':
            context['confirm_to_upload']=self.storage.extra_data.get('confirm_to_upload', None)
        return context
    