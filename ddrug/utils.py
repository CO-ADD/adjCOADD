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

from apputil.utils import Filterbase, get_filewithpath
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from adjcoadd.constants import *
from django.conf import settings
# ======================================Util Func. (To SVG)=====================================================#


def molecule_to_svg(mol, file_name, width=500, height=500):
    """Save substance structure as SVG"""
   
    file_path=get_filewithpath(file_name=file_name)
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
    for filename in os.listdir(path):
        file_path=os.path.join(path, filename)
        try:
            os.unlink(file_path)
            print("removed!")
        except Exception as err:
            print(err)


# def smiles_substructure_query(substructure):
#     query = Compound.objects.filter(molecule__hassubstruct=substructure)
#     for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
#         print(cmpd.name, cmpd.smiles)

# molecular choices:
mol_choices=[(),(),()]


class Drug_filter(Filterbase):
    drug_name = django_filters.CharFilter(lookup_expr='icontains')
    smol=django_filters.CharFilter(lookup_expr='contains', ) #method='substructure_filter',)
    class Meta:
        model=Drug
        fields=['drug_name', 'smol']



class Vitekcard_filter(Filterbase):
    card_barcode = django_filters.CharFilter(lookup_expr='icontains')
    class Meta:
        model=VITEK_Card
        fields=['card_barcode']