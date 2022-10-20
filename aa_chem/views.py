import re
import csv
import os
from os.path import exists
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.core.paginator import Paginator
from django.core import management
from django.shortcuts import HttpResponse, render, redirect
from django_rdkit.models import * 
from aa_chem.models import Drugbank,Taxonomy, Organisms
from app.models import Dictionaries
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
from .forms import CreateNewOrgForm
from model_utils import Choices
from django.forms import modelform_factory
# ======================================Util Func. (To SVG)=====================================================#

def molecule_to_svg(mol, file_name, width=500, height=500):
    """Save substance structure as Png"""

    # Define full path name
    full_path = f"static/images/{file_name}.svg"

    # Render high resolution molecule
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Export to png
    cairosvg.svg2svg(bytestring=drawer.GetDrawingText().encode(), write_to=full_path)

# =======================================Taxo View=============================================================================#

# Create your views here.
@login_required
def home(req):
    for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)
    # search function
    if req.method=='POST':
        search =req.POST.get('search')
        field=req.POST.get('field')
        if field=='Organism_Name':
            result=Taxonomy.objects.filter(Organism_Name__contains=search)
        # elif field=='status':
        #     result=Taxonomy.objects.filter(status__contains=search)
        # else:
        #     result=Taxonomy.objects.filter(drug_id__contains=search)
        #     print(result)
    else:
        result=Taxonomy.objects.all()


    # objects_all=Taxonomy.objects.all()
    objects_all=result
    p=Paginator(objects_all, 24)
    page_number = req.GET.get('page')
    page_obj=p.get_page(page_number)
    for object_ in page_obj:
       
            # if exists(f"static/images/{object_.id}.png"):
            #     print("file exists")
            #     return context
            # else:

        # if object_.id%2==0:
        m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
        # else:
            # m=Chem.MolFromSmiles('CCCC[C@@H]1NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](N)Cc2ccccc2)C(C)C)CCC(=O)NCCCC[C@@H](C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H](C(=O)C(N)=O)[C@@H](C)CC)[C@@H](C)CC)NC(=O)[C@H](C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](C)NC1=O')
        molecule_to_svg(m, object_.Organism_Name)
            
 

   
    context={
        'page_obj':page_obj,
        'chose':Taxonomy.Choice_Dictionaries   
       
    }
  
    return render(req, 'aa_chem/chem.html', context)

class TaxoListView(ListView):
    model=Taxonomy
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/taxoListview.html'

    def get_context_data(self, **kwargs):
        for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)

        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        objects_all=[object_ for object_ in context["objects"]]
        p=Paginator(objects_all, 24)
        
        for object_ in p.get_page(self.request.GET.get('page')):
       
            # if exists(f"static/images/{object_.id}.png"):
            #     print("file exists")
            #     return context
            # else:

            # if object_.id%2==0:
            m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
            # else:
                # m=Chem.MolFromSmiles('CCCC[C@@H]1NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](N)Cc2ccccc2)C(C)C)CCC(=O)NCCCC[C@@H](C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H](C(=O)C(N)=O)[C@@H](C)CC)[C@@H](C)CC)NC(=O)[C@H](C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](C)NC1=O')
            molecule_to_svg(m, object_.Organism_Name)
            
        return context



class TaxoCreateView(CreateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoCreate.html'
    success_url = reverse_lazy('compounds')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
    
        context["objects"]=self.model.objects.all()
        return context

class TaxoUpdateView(UpdateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoUpdate.html'
    success_url = reverse_lazy('compounds')

# ===============================================================OrgView==============================================#
# class OrgCreateView(CreateView):
#     model=Organisms
#     form_class=CreateNewOrgForm
#     template_name = 'aa_chem/orgCreate2.html'
#     success_url = reverse_lazy('compounds')
OrgForm=modelform_factory(Organisms, exclude=["Strain_Type"])
def newOrgnisms(req):
    org_strain=Dictionaries.objects.filter(Dictionary_ID='Units_Concentration')
    strains=str(org_strain[0])
    strain_test=strains.split(" ")[0].split(',')
    print(strain_test)
    a=strain_test[0]
    b=strain_test[0]
    c=strain_test[0]

    if req.method=='POST':
        form=OrgForm(req.POST)
        strain=req.POST.get('strain')
        try:
            print('iam trying')
            if form.is_valid():
                instance=form.save()
                # instance.save()
                instance.Strain_Type=str(strain)
                instance.save(update_fields=['Strain_Type'])
                print('saved!')
                return redirect("/")
        except Exception as err:
            print(err)
    else:
        form=OrgForm()
    return render(req, 'aa_chem/orgCreate2.html', {'form':form, 'a':a, 'b':b, 'c':c})


class OrgListView(ListView):
    model=Organisms
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/orgList.html'
    
    def get_context_data(self, **kwargs):
        for filename in os.listdir("static/images"):
                file_path=os.path.join("static/images", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)

        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        objects_all=[object_ for object_ in context["objects"]]
        p=Paginator(objects_all, 24)
        
        for object_ in p.get_page(self.request.GET.get('page')):
       
            # if exists(f"static/images/{object_.id}.png"):
            #     print("file exists")
            #     return context
            # else:

            if object_.id%2==0:
                m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
            else:
                m=Chem.MolFromSmiles('CCCC[C@@H]1NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](N)Cc2ccccc2)C(C)C)CCC(=O)NCCCC[C@@H](C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H](C(=O)C(N)=O)[C@@H](C)CC)[C@@H](C)CC)NC(=O)[C@H](C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](C)NC1=O')
            molecule_to_svg(m, object_.id)
            
        return context

class OrgTableView(OrgListView):
  
    template_name = 'aa_chem/orgTable.html'

   



# ==============================Function View===========================================================#


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
   