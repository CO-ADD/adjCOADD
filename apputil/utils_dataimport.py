import magic
import os
import pandas as pd
from pathlib import Path


from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.storage import FileSystemStorage
from django.template.defaultfilters import filesizeformat
from django.utils.deconstruct import deconstructible
from django.views import View
from django.shortcuts import render

from .utils import instance_dict, Validation_Log, SuperUserRequiredMixin


# -----------------------Start Utility Functions-----------------------------------

# ==============Uploading File Validators==============================

@deconstructible
class FileValidator(object):
    error_messages = {
     'max_size': ("Ensure this file size is not greater than %(max_size)s."
                  " Your file size is %(size)s."),
     'min_size': ("Ensure this file size is not less than %(min_size)s. "
                  "Your file size is %(size)s."),
     'content_type': "Files of type %(content_type)s are not supported.",
    }

    def __init__(self, max_size=None, min_size=None, content_types=()):
        self.max_size = max_size
        self.min_size = min_size
        self.content_types = content_types

    def __call__(self, data):
        # if self.max_size is not None and data.size > self.max_size:
        #     params = {
        #         'max_size': filesizeformat(self.max_size), 
        #         'size': filesizeformat(data.size),
        #     }
        #     raise ValidationError(self.error_messages['max_size'],
        #                            'max_size', params)

        # if self.min_size is not None and data.size < self.min_size:
        #     params = {
        #         'min_size': filesizeformat(self.min_size),
        #         'size': filesizeformat(data.size)
        #     }
        #     raise ValidationError(self.error_messages['min_size'], 
        #                            'min_size', params)

        if self.content_types:
            content_type = magic.from_buffer(data.read(), mime=True)
            data.seek(0)

            if content_type not in self.content_types:
                params = { 'content_type': content_type }
                raise ValidationError(self.error_messages['content_type'],
                                   'content_type', params)

    def __eq__(self, other):
        return (
            isinstance(other, FileValidator) and
            self.max_size == other.max_size and
            self.min_size == other.min_size and
            self.content_types == other.content_types
        )

# ==========================File Process=====================================

# set filefield Validator
validate_file = FileValidator(#max_size=1024 * 100, 
                             content_types=('text/csv', 'application/pdf','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
# create array for files if infected
# infected_files = []
# setup unix socket to scan stream
# cd = clamd.ClamdUnixSocket()




class FileUploadForm(SuperUserRequiredMixin, forms.Form):
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True, 'webkitdirectory':True}), validators=[validate_file])
    

class Importhandler(SuperUserRequiredMixin, View):
    
    form_class=FileUploadForm
    file_url=[]
    data_list=[]
    data_model='default'
    success_url="default"
    template_name='default'
    # log_process='default'
    vLog = Validation_Log("")
    # validatefile_name=["|"]
    validate_result={}
    file_report={}
    dirname=settings.MEDIA_ROOT

    #---------------------------------------------------------------------------------- 
    # Use to delete uploaded files
    def delete_file(self, file_name):
       
        file_full_path=os.path.join(self.dirname, file_name)
        print(file_full_path)
        try:
            os.unlink(file_full_path)
            print("removed!")
        except Exception as err:
            print(err)

    #---------------------------------------------------------------------------------- 
    # use to validate records
    def validates (self, newentry_dict, app_model, vlog, report_result, report_filelog, save, **kwargs):
       
        if any(newentry_dict.values()):
            for key in newentry_dict:
                for e in newentry_dict[key]:
                    djCard=app_model.check_from_dict(e, vlog)
                    if djCard.validStatus:
                        if save:
                            try:
                                djCard.save(**kwargs)
                                e['validStatus']=True
                            except Exception as err:
                                self.validate_result[key].append(f"catch Exception Card {err}")
                        else:
                            e['validStatus']=False
                    if report_result.get(key, None):
                        report_result[key].append(str(app_model._meta.db_table)+'_'+str(djCard.validStatus))
                    else:
                        report_result[key]=[str(app_model._meta.db_table)+'_'+str(djCard.validStatus)]
                    if report_filelog.get(key, None):
                        report_filelog[key].append(str(vlog.show()))
                    else:
                        report_filelog[key]=[str(vlog.show())]
                if save:
                    newentry_dict[key][:]=[e for e in newentry_dict[key] if e['validStatus']==False]
                    
    #---------------------------------------------------------------------------------- 
    def get(self, request):
        form = self.form_class
        template=self.template_name
        chars_lookup=str(request.user)
        # filesinUploads_list=[file for file in os.listdir(self.dirname) if os.path.isfile(os.path.join(self.dirname, file)) and str(request.user)+"_" in file] 
               
        return render(request, 'ddrug/importhandler_vitek.html', { 'form': form,})




