from django.http import JsonResponse
from .models import Organism, Taxonomy

#--Ajax search funcion--
## Search Organism Name in Taxonomy
def search_organism(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        searchInput=req.POST.get('inputtext')
        qs=Taxonomy.objects.filter(organism_name__icontains=searchInput)
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

## Search Organism ID in Organism
def search_organism_id(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        searchInput=req.POST.get('inputtext')
        qs=Organism.objects.filter(organism_id__icontains=searchInput)
        if len(qs)>0 and len(searchInput)>0:
            data=[]
            for i in qs:
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

