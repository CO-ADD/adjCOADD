"""
File-upload...
"""
import os
import magic

from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.storage import FileSystemStorage
from django.template.defaultfilters import filesizeformat
from django.utils.deconstruct import deconstructible
from django.views import View
from django.shortcuts import render

from .views_base import SuperUserRequiredMixin


# --files upload--
## create user folder
def file_location(req):
    location=settings.MEDIA_ROOT+'/'+str(req.user)
    return location

## Override filename in FileStorage
class OverwriteStorage(FileSystemStorage):
    
    def get_available_name(self, name, max_length=None):
        """
        Returns a filename that's free on the target storage system, and
        available for new content to be written to.
        """
        # If the filename already exists, remove it as if it was a true file system
    
        if self.exists(name):
            os.remove(os.path.join(self.location, name))
        return name

## validate files size, type, scan virus...
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

# 

# set filefield Validator
validate_file = FileValidator(#max_size=1024 * 100, 
                             content_types=('text/csv', 'application/pdf','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'))
# create array for files if infected
# infected_files = []
# setup unix socket to scan stream
# cd = clamd.ClamdUnixSocket()
##
class MultiFileUploadForm(SuperUserRequiredMixin, forms.Form):
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True, 'webkitdirectory':True}), validators=[validate_file])

##
class FileUploadForm(SuperUserRequiredMixin, forms.Form):
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(),  validators=[validate_file])

##
class Importhandler(SuperUserRequiredMixin, View):  
    """
    upload, parse and import data from pdf
    """  
    form_class=MultiFileUploadForm
    file_url=[]
    data_model='default'
    validate_result={}
    file_report={}
    process_name=None
    # Use to delete uploaded files
    def delete_file(self, file_name):
        location=file_location(self.request)
        file_full_path=os.path.join(location, file_name)
        print(file_full_path)
        try:
            os.unlink(file_full_path)
            print("removed!")
        except Exception as err:
            raise Exception
    # use to validate records
    def validates (self, newentry_dict, app_model, vlog, report_result, report_filelog, save, **kwargs):
        if any(newentry_dict.values()):
            for key in newentry_dict:
                for e in newentry_dict[key]:
                    newentry=app_model.check_from_dict(e, vlog)
                    if newentry.VALID_STATUS:
                        if save:
                            try:
                                newentry.save(**kwargs)
                                e['VALID_STATUS']=True
                            except Exception as err:
                                self.validate_result[key].append(f"catch Exception Card {err}")
                        else:
                            e['VALID_STATUS']=False
                    if report_result.get(key, None):
                        report_result[key].append(str(app_model._meta.db_table)+'_'+str(newentry.VALID_STATUS))
                    else:
                        report_result[key]=[str(app_model._meta.db_table)+'_'+str(newentry.VALID_STATUS)]
                    if report_filelog.get(key, None):
                        report_filelog[key].append(str(vlog.log_to_UI()))
                    else:
                        report_filelog[key]=[str(vlog.log_to_UI())]
                    
    def get(self, request, process_name):
        self.process_name=process_name
        return render(request, self.template_name, { 'form': self.form_class, 'process_name':self.process_name})



