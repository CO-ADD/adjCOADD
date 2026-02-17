import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError
from dsample.models import COADD_Compound, Compound_Batch
from dsample.utils.compounds import Upload_COADD_Compound


class Command(BaseCommand):
    help = "Uploads COADD samples to Project"

    def add_arguments(self, parser):
        parser.add_argument("--csv_file")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        
    def handle(self, *args, **options):
        if options["csv_file"]:
            if os.path.isfile(options["csv_file"]):
                samples = pd.read_csv(options["csv_file"])
                
                for idx,row in samples.iterrows():
                    Upload_COADD_Compound(row, upload=options["upload"], overwrite=options["overwrite"], appuser=None)
                

