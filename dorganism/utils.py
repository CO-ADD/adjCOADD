import os
import django_filters

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from .models import Organism, Taxonomy
from apputil.models import Dictionary



# ===================================Multiple choices__Dictionary queryset convert to choice Tuples ==#

def querysetToChoiseList_Dictionary(model_name, field_name):

    options=model_name.objects.filter(dict_class=field_name).values('dict_value', 'dict_desc')
    if options:
       
        choices_test=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]+ " | "+a[1]) for a in choices_test)
       
    else:
        choices=(('--', 'empty'),)
    return choices



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
    org_class=django_filters.ChoiceFilter(field_name='org_class', choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['organism_class']))


    class Meta:
        model=Taxonomy
        fields=['organism_name']

    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)


class Organismfilter(Filterbase):
    # pass
    organism_name = django_filters.CharFilter(field_name='organism_name__organism_name', lookup_expr='icontains')
    organism_class=django_filters.ChoiceFilter(field_name='organism_name__org_class__dict_value', choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['organism_class']))
    strain_type=django_filters.MultipleChoiceFilter(method='my_custom_filter', choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_type']))
    class Meta:
        model=Organism
        fields=['organism_id', 'strain_code', 'strain_ids',  'strain_notes', 'strain_type', 'mta_document', ]
       
    def my_custom_filter(self, queryset, name, value):
        return queryset.filter(strain_type__overlap=value)


class Taxonomyfilter(Filterbase):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    lineage = django_filters.MultipleChoiceFilter(choices="")
    class Meta:
        model=Taxonomy
        fields=['organism_name', 'code', 'org_class', 'tax_id', 'parent_tax_id', 'tax_rank', 'division', 'lineage']