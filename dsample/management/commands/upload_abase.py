import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log

from dsample.models import Project
from dsample.utils.compounds import Upload_COADD_Compound

#from oraABase.oraABase import openABase
from applib.mol.abase import get_ABase_RegView, upload_ABase_RegView, get_ABase_Batches, upload_ABase_Batches

from dcollab.models import Collab_Group, Collab_User, Collab_Membership, Organisation

from applib.project.import_project import get_CompoundSubmisssion_xlsx, parse_SampleInfo_Sheet, parse_ContactInfo_Sheet, Upload_Project_Collab


class Command(BaseCommand):
    help = "Uploads ABase compound batches"

    def add_arguments(self, parser):
        # parser.add_argument("--projectid")
        # parser.add_argument("--xlsx_file")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Upload_Project")

        #dfABase = get_ABase_RegView(test=0)
        #upload_Abase_RegView(dfABase, upload=options['upload'], overwrite=options['overwrite'])
        
        dfABase = get_ABase_Batches(test=0)
        print(f' [upload_abase] Projects : {dfABase['study_id'].unique()}')
        #upload_ABase_Batches(dfABase, upload=options['upload'], overwrite=options['overwrite'])