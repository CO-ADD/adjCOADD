import re
import csv
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.core.paginator import Paginator
from django.core import management
from django.shortcuts import HttpResponse, render, redirect
from django_rdkit.models import * 
from aa_chem.models import Drugbank,Taxonomy
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
import cairosvg
import py3Dmol
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.urls import reverse_lazy



class TaxoCreateView(CreateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoCreate.html'
    success_url = reverse_lazy('compounds')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context

class TaxoListView(ListView):
    model=Taxonomy
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/taxoList.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        
        print(context["objects"])
        
        return context



def molecule_to_png(mol, file_name, width=1000, height=1000):
    """Save substance structure as Png"""

    # Define full path name
    full_path = f"static/images/{file_name}.png"

    # Render high resolution molecule
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
    cairosvg.svg2png(bytestring=drawer.GetDrawingText().encode(), write_to=full_path)

# Create your views here.
@login_required
def home(req):

    # search function
    if req.method=='POST':
        search =req.POST.get('search')
        field=req.POST.get('field')
        if field=='drug_name':
            result=Drugbank.objects.filter(drug_name__contains=search)
        elif field=='status':
            result=Drugbank.objects.filter(status__contains=search)
        else:
            result=Drugbank.objects.filter(drug_id__contains=search)
            print(result)
    else:
        result=Drugbank.objects.all()
 
    # comp_all=[comp for comp in Drugbank.objects.annotate(amw=AMW('drug_mol'))]
    comp_all=[comp for comp in Drugbank.objects.annotate(smiles=MOL_TO_SMILES('drug_mol'))]

   ##!!will move this function to dataimport  
    for comp in comp_all:
        
        m=Chem.MolFromSmiles(comp.smiles)
        molecule_to_png(m, comp.id)
    

    paginator = Paginator(comp_all, 2) #3 elements per page.
    page_number = req.GET.get('page')
    comp_all = paginator.get_page(page_number)

   
    context={
        'mylist':comp_all,
       
       
    }
  
    return render(req, 'aa_chem/chem.html', context)




# @staff_member_required
# @user_passes_test(lambda u: u.is_staff)
@permission_required('app.change_groupfilter')
# @login_required
def importCSV(req):  

    if req.method=='POST':
        try:
            fdata=[req.POST.get('fdata'),]
            management.call_command('processCsvmultiple', fpath=fdata)
            print("done!")
            print(req.body)
        except Exception as err:
            print(err)
   
    return redirect(req.META['HTTP_REFERER'])

    
# @user_passes_test(lambda u: u.is_admin)
@permission_required('importdata')
def exportCSV(req):
    response=HttpResponse(content_type='text/csv')
   
    writer=csv.writer(response)
    writer.writerow(['S', 'O','C', "Cl", "NC", "NP", "T","D", "Di", "Lineaage"])#['id', 'drug_name','drug_mol']
    query=Taxonomy.objects.all()
    comp_list=[comp for comp in query]

    for comp in comp_list:
        writer.writerow(['S', 'O','C', "Cl", "NC", "NP", "T","D", "Di", "Lineaage"]) #[comp.id, comp.drug_name, comp.smiles]

    response['Content-Disposition']='attachment; filename="Taxo_export.csv"'
    return response
   