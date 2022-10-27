from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem.Draw import IPythonConsole
from IPython.display import SVG
import cairosvg
# import py3Dmol
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


#=================================================Clear IMGFolder===========================================================#

def clearIMGfolder():
    for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)

#===================== A Example =================================================================================#
# def home(req): 
#     clearIMGfolder()
#     # search function
#     if req.method=='POST':
#         search =req.POST.get('search')
#         field=req.POST.get('field')
#         if field=='Organism_Name':
#             result=Taxonomy.objects.filter(astatus=1, Organism_Name__contains=search)
#     else:
#         result=Taxonomy.objects.filter(astatus=1)
        
#     objects_all=result
#     p=Paginator(objects_all, 24)
#     page_number = req.GET.get('page')
#     page_obj=p.get_page(page_number)
#     for object_ in page_obj:
#         m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
#         molecule_to_svg(m, object_.Organism_Name)
    
#     context={
#         'page_obj':page_obj,
#         'chose':Taxonomy.Choice_Dictionaries   
       
#     }
  
#     return render(req, 'aa_chem/chem.html', context)

# ===================================Dictionary query convert to choice Tuples========================================================================#


def querysetToChoiseList_Dictionaries(model_name, field_name):
    options=model_name.objects.filter(Dictionary_Class=field_name).values('Dict_Value', 'Dict_Desc')
    if options:

        choices=[tuple(d.values()) for d in options]
    else:
        choices=(('--', 'empty'),)
    return choices