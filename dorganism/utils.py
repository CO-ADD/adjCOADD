import os
import django_filters

from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from .models import Organism, Taxonomy, Organism_Batch
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


def search_organism_id(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        searchInput=req.POST.get('inputtext')
        print(searchInput)
        qs=Organism.objects.filter(organism_id__istartswith=searchInput)
        print(qs)
        if len(qs)>0 and len(searchInput)>0:
            data=[]
            for i in qs:
                print(i.organism_id)
                item={
                    'name':i.organism_id,
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
   
    def multichoices_filter(self, queryset, name, value):
        lookup='__'.join([name, 'overlap'])
        return queryset.filter(**{lookup: value})


class Organismfilter(Filterbase):
    organism_id=django_filters.CharFilter(field_name='organism_id', lookup_expr='icontains')
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
       



class Taxonomyfilter(Filterbase):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    lineage = django_filters.MultipleChoiceFilter( choices= "")
    # django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['lineage'], ' | '))
    division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    class Meta:
        model=Taxonomy
        fields=['organism_name', 'code', 'org_class', 'tax_id', 'parent_tax_id', 'tax_rank', 'division', 'lineage']




class Batchfilter(Filterbase):
    Stock_Date=django_filters.IsoDateTimeFilter(field_name='stock_date')
    class Meta:
        model=Organism_Batch
        fields= ["supplier","supplier_code","supplier_po", "stock_date",  "biologist"]



# -------------------editable tables utility function--------------------------------
# @user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
# @csrf_protect
# def detailChangeOrganism(req):
#     kwargs={}
#     kwargs['user']=req.user 
#     id=req.POST.get('id', '')
#     object_=get_object_or_404(Organism, organism_id=id)
#     value=req.POST.get('value','')
#     type_value=req.POST.get('type', '')

#     if type_value=='strain_type':
#         try:
#             value=value.split(",")
#             object_.strain_type=[i for i in value]
#             object_.save(**kwargs)
#         except Exception as err:
#              print("something wroing")
    
#     else:
#         try:
#             fields={type_value: value}
#             print(fields)
#             Organism.objects.filter(pk=id).update(**fields)
#             object_=get_object_or_404(Organism, organism_id=id)
#             object_.save(**kwargs)
#         except Exception as err:
#             print(err)
   
#     return JsonResponse({"success": "updated!"})