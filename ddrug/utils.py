import os
from pathlib import Path
# os.environ['path']+=r';C:\Program Files\UniConvertor-2.0rc5\dlls'
import django_filters
from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg
# import py3Dmol

from dorganism.utils import Filterbase
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from adjcoadd.constants import *
from django.conf import settings
# ======================================Util Func. (To SVG)=====================================================#
def molecule_to_svg(mol, file_name, width=500, height=500):
    """Save substance structure as SVG"""
   
    # Define full path name
    if settings.DEVELOPMENT:
        file_path = f"static/images/{file_name}.svg"
    # print(f'path1: {file_path1}')
    else:
        Base_dir = Path(__file__).resolve().parent.parent.parent
        FILES_DIR=os.path.abspath(os.path.join(Base_dir, 'static/images'))
        file_path=os.path.join(FILES_DIR, f"{file_name}.svg") 

    # Render high resolution molecule
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
 
    cairosvg.svg2svg(bytestring=drawer.GetDrawingText().encode(), write_to=file_path)
  

#=================================================Clear IMGFolder===========================================================#

def clearIMGfolder():
    if settings.DEVELOPMENT:
        path='static/images'
    else:
        Base_dir = Path(__file__).resolve().parent.parent.parent
        path=os.path.abspath(os.path.join(Base_dir, 'static/images'))
    for filename in os.listdir(path): # os.listdir("static/images/"):
        file_path=os.path.join(path, filename)
        try:
            os.unlink(file_path)
            print("removed!")
        except Exception as err:
            print(err)


class Drug_filter(Filterbase):
    drug_name = django_filters.CharFilter(lookup_expr='icontains')
    # lineage = django_filters.MultipleChoiceFilter( choices= "")
    # django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['lineage'], ' | '))
    # division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    class Meta:
        model=Drug
        fields=['drug_name']



class Vitekcard_filter(Filterbase):
    card_barcode = django_filters.CharFilter(lookup_expr='icontains')
    # lineage = django_filters.MultipleChoiceFilter( choices= "")
    # django_filters.MultipleChoiceFilter(method='multichoices_filter', choices=get_DictonaryChoices_byDictClass(Dictionary, Organism.Choice_Dictionary['lineage'], ' | '))
    # division= django_filters.ModelChoiceFilter(queryset=Dictionary.objects.filter(dict_class=Taxonomy.Choice_Dictionary['division']))
    class Meta:
        model=VITEK_Card
        fields=['card_barcode']