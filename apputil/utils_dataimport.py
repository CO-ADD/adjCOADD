import magic
import magic
import os
import pandas as pd
from pathlib import Path


from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.storage import FileSystemStorage
from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.core.files.storage import FileSystemStorage
from django.db import transaction
from django.template.defaultfilters import filesizeformat
from django.utils.deconstruct import deconstructible
from django.views.generic.edit import FormView
from django.template.defaultfilters import filesizeformat
from django.utils.deconstruct import deconstructible
from django.views.generic.edit import FormView

from dorganism.models import Taxonomy, Organism, Organism_Batch, Organism_Culture
from apputil.models import Dictionary
from apputil.utils import Validation_Log, instance_dict
from ddrug.models import VITEK_Card, VITEK_ID, VITEK_AST
from .utils import instance_dict, Validation_Log, SuperUserRequiredMixin
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
    
    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True,}), validators=[validate_file])
    

class Importhandler(SuperUserRequiredMixin, FormView):
    
    form_class=FileUploadForm
    file_url=[]
    data_list=[]
    data_model='default'
    success_url="default"
    template_name='default'
    log_process='default'
    vLog = Validation_Log(log_process)
    table_name=[","]
    validate_result={}
    file_report=[","]
    
    def delete_file(self, file_path):
        file_name=file_path.split("/")[2]
        print(file_name)
        file_full_path=os.path.join(settings.MEDIA_ROOT, file_name)
        print(file_full_path)
        try:
            os.unlink(file_full_path)
            print("removed!")
        except Exception as err:
            print(err)
 




