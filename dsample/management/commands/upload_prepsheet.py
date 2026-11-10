import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from apputil.models import ApplicationUser
from dsample.models import Project
from dsample.process_views import Upload_StockPrep_Process

from dsample.utils.compounds import Upload_COADD_Compound
from dcollab.models import Collab_Group, Collab_User, Collab_Membership, Organisation

from applib.plate.stockprep import get_StockPrep_xlsx

class Command(BaseCommand):
    help = "Uploads PrepSheet (Barcodes) to Project"

    def add_arguments(self, parser):
        parser.add_argument("--projectid")
        parser.add_argument("--xlsx")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Upload_PrepSheet")
        
        cmduser = ApplicationUser.get('J.Zuegg')
        if options["xlsx"] and options["projectid"]:
            # - Project Init ------------------------------------------------------
            djPrj = Project.get(options["projectid"])
            if djPrj:
                if os.path.isfile(options["xlsx"]):
                    dirname, filename = os.path.split(options["xlsx"])


                    valLog=Upload_StockPrep_Process(None, dirname, [filename], ProjectID=options["projectid"], 
                                                upload=options["upload"], overwrite=options["overwrite"],
                                                appuser=cmduser)
                else:
                    print(f"ERROR : PrepSheet {options["xlsx_file"]} does not exists")        
            else:
                print(f"ERROR : Project {options["projectid"]} does not exists")        
                    
        print("---------------------------------------------")
        valLog.select_unique()        
        valLog.show()
        print("---------------------------------------------")        
