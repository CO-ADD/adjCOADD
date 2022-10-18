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
            csv_reader = csv.reader(csv_file, delimiter=';')
            next(csv_reader) #advance past the header
            
            #Taxonomy.objects.all().delete()
            
            for row in csv_reader:
                line=row[9].split(",")
                print(row)
                obj, created=Taxonomy.objects.get_or_create(Organism_Name=row[0], Other_Names=row[1], Code=row[2], Class=row[3], Tax_ID=row[4],Parent_Tax_ID=row[5], Tax_Rank=row[6], Division=row[7], Division_CODE=row[8], Lineage=line)
