import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log

from dsample.models import Project
from dsample.utils.compounds import Upload_COADD_Compound
from dcollab.models import Collab_Group, Collab_User, Collab_Membership, Organisation

from applib.project.import_project import get_CompoundSubmisssion_xlsx, parse_SampleInfo_Sheet, parse_ContactInfo_Sheet, Upload_Project_Collab

class Command(BaseCommand):
    help = "Uploads COADD samples to Project"

    def add_arguments(self, parser):
        parser.add_argument("--projectid")
        parser.add_argument("--xlsx_file")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Upload_Project")
        
        if options["xlsx_file"] and options["projectid"]:
        
            if os.path.isfile(options["xlsx_file"]):
                _dictSheets = get_CompoundSubmisssion_xlsx(options["xlsx_file"],valLog=valLog)
                
                
                # - Project Init ------------------------------------------------------
                djPrj = None
                if options["projectid"]:
                    djPrj = Project.get(options["projectid"])
                    if djPrj is None:
                        djPrj = Project()
                        djPrj.project_id = options["projectid"]
                        valLog.add_warning(f"New Project",f"{options["projectid"]}")
                    else:
                        valLog.add_info(f"Existing Project",f"{djPrj}")
                else:
                    djPrj = Project()
                    valLog.add_warning(f"New Project",f"{options["projectid"]}")

                # - Contacts and Project ----------------------------------------------
                if _dictSheets['Contacts'] is not None:
                    _Contacts,_PrjTitle = parse_ContactInfo_Sheet(_dictSheets['Contacts'],valLog=valLog)

                    djPrj.project_name = _PrjTitle
                    Upload_Project_Collab(djPrj, _Contacts, upload=options["upload"], overwrite=options["overwrite"], valLog=valLog)

                #- Samples ----------------------------------------------
                if _dictSheets['Samples'] is not None:
                    _Samples = parse_SampleInfo_Sheet(_dictSheets['Samples'],valLog=valLog)
                    for _smp in _Samples:
                        _smp['project_id'] = str(djPrj)
                        Upload_COADD_Compound(_smp, upload=options["upload"], overwrite=options["overwrite"], valLog=valLog)

                 
                    
                    #print(_Contacts)
                    
        print("---------------------------------------------")
        valLog.select_unique()        
        valLog.show()
        print("---------------------------------------------")        
