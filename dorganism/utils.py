import os
import django_filters
import pandas as pd

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from .models import Organism, Taxonomy
from apputil.models import Dictionary
from apputil.utils import get_DictonaryChoices_byDictClass

#=====================================search_organism===============================================================

def search_organism(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        searchInput=req.POST.get('inputtext')
        print(searchInput)
        qs=Taxonomy.objects.filter(organism_name__istartswith=searchInput)
        print(qs)
        if len(qs)>0 and len(searchInput)>0:
            data=[]
            for i in qs:
                if i.org_class:
                    orgClass=i.org_class.dict_value
                else:
                    orgClass='noClass by Import or ...'
                
                item={
                    'name':i.organism_name,
                    'class': orgClass,
                }
                data.append(item)
            res=data
        else:
            res='No organism found...'
        
        return JsonResponse({'data':res})
    return JsonResponse({})



#==================================Filters======================================
class Filterbase(django_filters.FilterSet):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    org_class=django_filters.ChoiceFilter(field_name='org_class', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['organism_class'], ' | '))


    class Meta:
        model=Taxonomy
        fields=['organism_name']

    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)


class Organismfilter(Filterbase):
    organism_name = django_filters.CharFilter(field_name='organism_name__organism_name', lookup_expr='icontains')
    organism_class=django_filters.ChoiceFilter(field_name='organism_name__org_class__dict_value', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['organism_class']))
    strain_type=django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['strain_type'], ' | '))
    strain_panel=django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['strain_panel'], ' | '))
    risk_group=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['risk_group']))
    mta_status=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['mta_status']))
    oxygen_pref=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['oxygen_pref']))
    pathogen_group=django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Organism.Choice_Dictionary['pathogen_group']))
    class Meta:
        model=Organism
        fields=['organism_id', 'strain_code', 'strain_ids', 'strain_type', 'mta_document', 'strain_panel', 'risk_group', 'mta_status', 'oxygen_pref', 'pathogen_group', ]
       
    def multichoices_filter(self, queryset, name, value):
        lookup='__'.join([name, 'overlap'])
        return queryset.filter(**{lookup: value})



class Taxonomyfilter(Filterbase):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    lineage = django_filters.MultipleChoiceFilter(choices="")
    division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    class Meta:
        model=Taxonomy
        fields=['organism_name', 'code', 'org_class', 'tax_id', 'parent_tax_id', 'tax_rank', 'division', 'lineage']

from django.db import transaction

def has_special_char(text: str) -> bool:
    return any(c for c in text if not c.isalnum() and not c.isspace())

@transaction.atomic
def import_excel(file_path):
    print('importing....')
    Taxonomy_object_list=[]
    excel_file=file_path
    print(excel_file)
    exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
    print(type(exmpexceldata))
    dbframe=exmpexceldata
    for dbframe in dbframe.itertuples():
        try:
            print(str(dbframe.ORGANISM_CLASS))
            if has_special_char(str(dbframe.ORGANISM_NAME)):
                print("speical")
                raise Exception('datatype')
            class_fkey=Dictionary.objects.filter(dict_value=dbframe.ORGANISM_CLASS)
            if class_fkey:
                class_fkey=class_fkey[0]
            else:
                class_fkey=None
            print(class_fkey)
            division_fkey=Dictionary.objects.filter(dict_value=dbframe.DIVISION)
            if division_fkey:
                division_fkey=division_fkey[0]
            else:
                division_fkey=None
            linea=str(dbframe.LINEAGE).split(";")
            print(division_fkey)
            try:
                Taxonomy_object_list.append(Taxonomy(organism_name=dbframe.ORGANISM_NAME, other_names=dbframe.ORGANISM_NAME_OTHER, code=dbframe.ORGANISM_CODE, 
                    org_class=class_fkey, tax_id=dbframe.TAX_ID, parent_tax_id=dbframe.PARENT_TAX_ID, 
                    tax_rank=dbframe.TAX_RANK, division=division_fkey, lineage=linea, 
                    ))
            except Exception as err:
                print(err)
                # obj.save()
            
        except Exception as err:
            print(err)
            return err
    return Taxonomy_object_list
    # return 'something wrong'