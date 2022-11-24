import os
import django_filters
from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg
# import py3Dmol
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

from .models import Organism, Taxonomy
from apputil.models import Dictionary


# ======================================Util Func. (To SVG)=====================================================#
def molecule_to_svg(mol, file_name, width=500, height=500):
    """Save substance structure as SVG"""

    # Define full path name
    file_path = f"static/images/{file_name}.svg"

    # Render high resolution molecule
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
    cairosvg.svg2svg(bytestring=drawer.GetDrawingText().encode(), write_to=file_path)


#=================================================Clear IMGFolder===========================================================#

def clearIMGfolder():
    for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)


# ===================================Dictionary query convert to choice Tuples========================================================================#


def querysetToChoiseList_Dictionary(model_name, field_name):

    options=model_name.objects.filter(dict_class=field_name).values('dict_value', 'dict_desc')
    if options:
       
        choices_test=tuple([tuple(d.values()) for d in options])
        choices=tuple((a[0], a[0]+' ( '+ a[1]+' )') for a in choices_test)
       
    else:
        choices=(('--', 'empty'),)
    return choices



#=====================================Searchbar_01===============================================================

def searchbar_01(req):
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


#==================================searchbar_02======================================
# def searchbar_02(req, model, model_field):
#     variable_column=model_field
#     search_type='contains'
#     filter_search=variable_column+'__'+search_type
#     # gt=None
#     if req.method=='POST':
#         qs =req.POST.get('search')
#         # gt=qs.strip()
#         field=req.POST.get('field')
#         if field==model_field:
#             result=model.objects.filter(astatus__gte=0, **{filter_search:qs})
        
#     else:
#         result=model.objects.filter(astatus__gte=0)
    
#     return result




#==================================searchbar_02======================================
class MySearchbar02(django_filters.FilterSet):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    class Meta:
        model=Taxonomy
        fields=['organism_name']

    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)


class MySearchbar03(MySearchbar02):
    # pass
    organism_name = django_filters.CharFilter(field_name='organism_name__organism_name', lookup_expr='icontains')
    organism_class=django_filters.ChoiceFilter(field_name='organism_name__org_class__dict_value', choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['organism_class']))
    strain_type=django_filters.MultipleChoiceFilter(method='my_custom_filter', choices=querysetToChoiseList_Dictionary(Dictionary, Organism.Choice_Dictionary['strain_type']))
    class Meta:
        model=Organism
        fields=['organism_id', 'strain_code', 'strain_id',  'strain_notes', 'strain_type', 'mta_document', ]
       
    def my_custom_filter(self, queryset, name, value):
        return queryset.filter(strain_type__overlap=value)


class MySearchbar04(MySearchbar02):
    organism_name = django_filters.CharFilter(lookup_expr='icontains')
    lineage = django_filters.MultipleChoiceFilter(choices="")
    class Meta:
        model=Taxonomy
        fields=['organism_name', 'code', 'org_class', 'tax_id', 'parent_tax_id', 'tax_rank', 'division', 'lineage']