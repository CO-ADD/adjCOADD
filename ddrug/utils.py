import os
from pathlib import Path
import django_filters
from django_rdkit.models import *
from django_rdkit.config import config
from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg

from django.conf import settings

from apputil.utils import Filterbase, get_filewithpath
from .models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID
from adjcoadd.constants import *
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


class Vitekast_filter(Filterbase):
    Drug_Name = django_filters.CharFilter(field_name='drug_id__drug_name', lookup_expr='icontains')
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.filters['Drug_Name'].label='Drug Name'
    class Meta:
        model=VITEK_AST
        fields=['Drug_Name']



# Similarity Query Function
# config.tanimoto_threshold =0.4 # similarity_threshold_int/100
def get_mfp2_neighbors(smiles):
    value = MORGANBV_FP(Value(smiles))
    print(config.tanimoto_threshold)
    # print(f'threshold {config.tanimoto_threshold}')
    queryset = Drug.objects.filter(mfp2__tanimoto=value)#(mfp2__tanimoto=value)
    # queryset = queryset.annotate(smiles=MOL_TO_SMILES('smol'))
    # queryset = queryset.annotate(smol=TANIMOTO_SML('mfp2', value))
    queryset = queryset.order_by(TANIMOTO_DIST('mfp2', value))
    # queryset = queryset.values_list('drug_name',  ) #'smiles','smol'
    return queryset