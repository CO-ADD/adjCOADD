import os
from pathlib import Path
from django_rdkit.models import *
from django_rdkit.config import config
from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg

from django.conf import settings

from ddrug.models import  Drug

#-----------------------------------------------------------------------------------
def smiles2mol(Smiles,verbose=0):
#-----------------------------------------------------------------------------------
    try:
        xmol = Chem.MolFromSmiles(Smiles)
    except:
        xmol = None
        if verbose:
            print(f"[Invalid SMILES] {Smiles} ")
    return(xmol)

# ----------------------------------------------------------------------------------------------------
# --Convert mol to Structure images--
## convert function
def molecule_to_svg(mol, file_name, width=500, height=500, path=settings.MOL_IMG_DIR):
    """Save substance structure as SVG"""
# ----------------------------------------------------------------------------------------------------
   
    file_svg=os.path.join(path, f"{file_name}.svg")  #get_filewithpath(file_name=file_name) 
    # Render high resolution molecule
   
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
 
    cairosvg.svg2svg(bytestring=drawer.GetDrawingText().encode(), write_to=file_svg)
  

# ----------------------------------------------------------------------------------------------------
## folder clean
def clearIMGfolder(path = settings.MOL_IMG_DIR,verbose=0):
# ----------------------------------------------------------------------------------------------------
    # if settings.DEVELOPMENT:
    #     path='static/images'
    # else:
    #     Base_dir = Path(__file__).resolve().parent.parent.parent
    #     path=os.path.abspath(os.path.join(Base_dir, 'static/images'))
    for filename in os.listdir(path):
        file_path=os.path.join(path, filename)
        try:
            os.unlink(file_path)
            print("removed!")
        except Exception as err:
            print(err)


# ----------------------------------------------------------------------------------------------------
# --Similarity Query Function--
# config.tanimoto_threshold =0.4 # similarity_threshold_int/100
def get_mfp2_neighbors(smiles):
# ----------------------------------------------------------------------------------------------------
    value = MORGANBV_FP(Value(smiles))
    print(config.tanimoto_threshold)
    # print(f'threshold {config.tanimoto_threshold}')
    queryset = Drug.objects.filter(mfp2__tanimoto=value)#(mfp2__tanimoto=value)
    # queryset = queryset.annotate(smiles=MOL_TO_SMILES('smol'))
    # queryset = queryset.annotate(smol=TANIMOTO_SML('mfp2', value))
    queryset = queryset.order_by(TANIMOTO_DIST('mfp2', value))
    # queryset = queryset.values_list('drug_name',  ) #'smiles','smol'
    return queryset

# def smiles_substructure_query(substructure):
#     query = Compound.objects.filter(molecule__hassubstruct=substructure)
#     for cmpd in query.annotate(smiles=MOL_TO_SMILES('molecule'))[:5]:
#         print(cmpd.name, cmpd.smiles)
