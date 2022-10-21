import os
import csv
from django.core.management.base import BaseCommand
from django.conf import settings
from aa_chem.models import Taxonomy, Organisms
from app.models import Dictionaries

class Command(BaseCommand):
   
    def handle(self, *args, **options):
        path=input("enter the path directory")
        with open(os.path.join(settings.BASE_DIR / path)) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader) #advance past the header
            
            
            
            for row in csv_reader:
                int_order=int('0'+row[5])
                dict_v=row[2].split(",")
                # print(row)
                obj, created=Dictionaries.objects.get_or_create(Dictionary_ID=row[0], Dictionary_Class=row[1], Dict_Value=dict_v, Dict_Desc =row[3], Dict_Value_Type =row[4], Dict_View_Order =int_order)

