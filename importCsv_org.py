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
                taxID=int('0'+row[22])
                screen_panel=row[26].split(';')
                try:
                    obj, created=Organisms.objects.get_or_create(Organism_ID=row[0], Organism_Name=row[1], Organism_Desc=row[2], Strain_ID=row[3], 
                                    Strain_Code=row[5], Strain_Desc=row[6], Strain_Notes=row[7], 
                                    Strain_Tissue=row[25], Strain_Type=row[4], Sequence=row[28], Sequence_Link=row[29], Geno_Type=row[33],
                                    Screen_Type=screen_panel, 
                                    Tax_ID =taxID,Risk_Group=row[9], Pathogen =row[10],Import_Permit =row[12],Biol_Approval =row[23],Special_Precaution =row[24],Lab_Restriction =row[27],MTA_Document =row[31],
                                    MTA_Status =row[32],Oxygen_Pref =row[13],Atmosphere_Pref =row[14],Nutrient_Pref =row[15],Biofilm_Pref =row[16], )
                
                except Exception as err:
                    print(err)



