from django.http import JsonResponse
from dorganism.models import Organism, Taxonomy


#-----------------------------------------------------------------------------------
def reformat_OrganismID(OrgID):
    """
    Reformat Organism ID from old GN_001 (3 digits) to new GN_0001 (3 digits) 
    """
#-----------------------------------------------------------------------------------
    xStr = OrgID.split("_")
    return(f"{xStr[0]}_{int(xStr[1]):04d}")

#-----------------------------------------------------------------------------------
def reformat_OrgBatchID(OrgBatchID):
    """
    Reformat OrganismBatch ID from old GN_001:02 (3 digits,':') to new GN_0001_02 (3 digits,'_') 
    """
#-----------------------------------------------------------------------------------
    xStr = OrgBatchID.split(":")
    #print(xStr)
    return(f"{reformat_OrganismID(xStr[0])}_{int(xStr[1]):02d}")

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

