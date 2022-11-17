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
from apputil.models import Dictionaries


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


def querysetToChoiseList_Dictionaries(model_name, field_name):
    options=model_name.objects.filter(Dictionary_Class=field_name).values('Dict_Value', 'Dict_Desc')
    if options:

        choices=tuple([tuple(d.values()) for d in options])
    else:
        choices=(('--', 'empty'),)
    return choices



#=====================================Searchbar_01===============================================================

def searchbar_01(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        searchInput=req.POST.get('inputtext')
        qs=Taxonomy.objects.filter(Organism_Name__istartswith=searchInput)
      
        if len(qs)>0 and len(searchInput)>0:
            data=[]
            for i in qs:
                if i.Class:
                    Class=i.Class.Dict_Value
                else:
                    Class='noClass by Import or ...'
                
                item={
                    'name':i.Organism_Name,
                    'class': Class,
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
    Organism_Name = django_filters.CharFilter(lookup_expr='icontains')
    class Meta:
        model=Taxonomy
        fields=['Organism_Name']

    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)


class MySearchbar03(MySearchbar02):
    # pass
    Organism_Name = django_filters.CharFilter(field_name='Organism_Name__Organism_Name', lookup_expr='icontains')
    Organism_Class=django_filters.ChoiceFilter(field_name='Organism_Name__Class__Dict_Value', choices=querysetToChoiseList_Dictionaries(Dictionaries, Organism.Choice_Dictionaries['Organism_Class']))
    Strain_Type=django_filters.MultipleChoiceFilter(method='my_custom_filter', choices=querysetToChoiseList_Dictionaries(Dictionaries, Organism.Choice_Dictionaries['Strain_Type']))
    class Meta:
        model=Organism
        fields=['Organism_ID', 'Strain_Code', 'Strain_ID',  'Strain_Notes', 'Strain_Type', 'MTA_Document', ]
       
    def my_custom_filter(self, queryset, name, value):
        return queryset.filter(Strain_Type__overlap=value)


class MySearchbar04(MySearchbar02):
    Organism_Name = django_filters.CharFilter(lookup_expr='icontains')
    Lineage = django_filters.MultipleChoiceFilter(choices="")
    class Meta:
        model=Taxonomy
        fields=['Organism_Name', 'Code', 'Class', 'Tax_ID', 'Parent_Tax_ID', 'Tax_Rank', 'Division', 'Lineage']