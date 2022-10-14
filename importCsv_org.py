import os
import csv
from django.core.management.base import BaseCommand
from django.conf import settings
from aa_chem.models import Taxonomy, Organisms
from app.models import *

class Command(BaseCommand):
   
    def handle(self, *args, **options):
        path=input("enter the path directory")
        with open(os.path.join(settings.BASE_DIR / path)) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader) #advance past the header
            
            
            
            for row in csv_reader:
                # print(row)
                obj, created=Organisms.objects.get_or_create(Organism_ID=row[0], Organism_Name=row[1], Organism_Desc=row[2], Strain_ID=row[3], 
                                Strain_Code=row[4], Strain_Desc=row[5], Strain_Notes=row[6], Strain_Tissue=row[7], Strain_Type=row[8], Sequence=row[9], sequence_Link=row[10], Geno_Type=row[11],
                                Screen_Type=row[12], Tax_ID =row[13],Risk_Group=row[14], Pathogen =row[15],Import_Permit =row[16],Biol_Approval =row[17],Special_Precaution =row[18],Lab_Restriction =row[19],MTA_Document =row[20],
                                MTA_Status =row[20],Oxygen_Pref =row[20],Atmosphere_Pref =row[20],Nutrient_Pref =row[20],Biofilm_Pref =row[20], )
