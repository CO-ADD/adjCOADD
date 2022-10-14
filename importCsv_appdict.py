import os
import csv
from django.core.management.base import BaseCommand
from django.conf import settings
from aa_chem.models import Taxonomy, Organisms
from app.models import ApplicationDictionary

class Command(BaseCommand):
   
    def handle(self, *args, **options):
        path=input("enter the path directory")
        with open(os.path.join(settings.BASE_DIR / path)) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            next(csv_reader) #advance past the header
            
            
            
            for row in csv_reader:
                int_order=int('0'+row[4])
                # print(row)
                obj, created=ApplicationDictionary.objects.get_or_create(dict_table=row[0], dict_field=row[1], dict_value=row[2], dict_value_type =row[3],dict_order =int_order, dict_desc=row[5])



