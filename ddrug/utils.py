import os
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
from .models import Drug
from adjcoadd.constants import *

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
    for filename in os.listdir("static/images/"):
                file_path=os.path.join("static/images/", filename)
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