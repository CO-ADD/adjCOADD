from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg
import py3Dmol
import os


# ======================================Util Func. (To SVG)=====================================================#
def molecule_to_svg(mol, file_name, width=500, height=500):
    """Save substance structure as SVG"""

    # Define full path name
    full_path = f"static/images/{file_name}.svg"

    # Render high resolution molecule
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
    cairosvg.svg2svg(bytestring=drawer.GetDrawingText().encode(), write_to=full_path)


#============================================================================================================#

def clearIMGfolder():
    for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)