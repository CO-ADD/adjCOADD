# import re
# import csv
# import os
import ast
# from os.path import exists
# from django.contrib.admin.views.decorators import staff_member_required
# from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.core.paginator import Paginator
# from django.core import management
from django.shortcuts import HttpResponse, render, redirect
# from django_rdkit.models import * 
from aa_chem.models import Mytest, Organisms, Taxonomy #Drugbank,,Organisms,
# import aa_chem.models
from app.models import Dictionaries, Mytest2
from app.utils import molecule_to_svg, clearIMGfolder
from rdkit import Chem

from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.urls import reverse_lazy
# # from .forms import CreateNewOrgForm
# from model_utils import Choices
# from django.forms import modelform_factory
from django.shortcuts import get_object_or_404
from django.http import JsonResponse




# # =======================================Taxo View=============================================================================#

# # Create your views here.
# @login_required
def home(req): 
    clearIMGfolder()
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
        
    objects_all=result
    p=Paginator(objects_all, 24)
    page_number = req.GET.get('page')
    page_obj=p.get_page(page_number)
    for object_ in page_obj:
       
            # if exists(f"static/images/{object_.id}.png"):
            #     print("file exists")
            #     return context
            # else:
        m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
        molecule_to_svg(m, object_.Organism_Name)
    
    context={
        'page_obj':page_obj,
        'chose':Taxonomy.Choice_Dictionaries   
       
    }
  
    return render(req, 'aa_chem/chem.html', context)




# class TaxoListView(ListView):
#     model=Taxonomy
#     paginate_by = 24
#     fields='__all__'
#     template_name = 'aa_chem/taxoListview.html'

#     def get_context_data(self, **kwargs):
#         for filename in os.listdir("static/images"):
#                 file_path=os.path.join("static/images", filename)
#                 try:
#                     os.unlink(file_path)
#                     print("removed!")
#                 except Exception as err:
#                     print(err)

#         context=super().get_context_data(**kwargs)
#         context["objects"]=self.model.objects.all()
#         objects_all=[object_ for object_ in context["objects"]]
#         p=Paginator(objects_all, 24)
        
#         for object_ in p.get_page(self.request.GET.get('page')):

#             m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
#             molecule_to_svg(m, object_.Organism_Name)
            
#         return context



class TaxoCreateView(CreateView):
    model=Taxonomy
    fields='__all__'
    template_name = 'aa_chem/taxoCreate.html'
    success_url = reverse_lazy('compounds')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
    
        context["objects"]=self.model.objects.all()
        return context

# class TaxoUpdateView(UpdateView):
#     model=Taxonomy
#     fields='__all__'
#     template_name = 'aa_chem/taxoUpdate.html'
#     success_url = reverse_lazy('compounds')

# # ===============================================================OrgView==============================================#
# # class OrgCreateView(CreateView):
# #     model=Organisms
# #     form_class=CreateNewOrgForm
# #     template_name = 'aa_chem/orgCreate2.html'
# #     success_url = reverse_lazy('compounds')

def home2(req):
    if req.method=='POST':
        #a) Add Note form
        tag=req.POST.get("note_tag")
            
        note_new= Mytest.objects.create(tag=tag)
        note_new.save()
        print('save')
        return redirect("/")
        
    return render(req, 'aa_chem/home2.html') 

# ============================Create new Organism======================================================#
# =============================1. Ajax Call search Taxo=============================================#

def searchTaxo(req):
    if req.headers.get('x-requested-with') == 'XMLHttpRequest':
        res=None
        taxo=req.POST.get('inputtext')
        # print(taxo)
        qs=Taxonomy.objects.filter(Organism_Name__icontains=taxo)
        if len(qs)>0 and len(taxo)>0:
            data=[]
            for i in qs:
                item={
                    'name':i.Organism_Name,
                }
                data.append(item)
            res=data
        else:
            res='No games found...'
        
        return JsonResponse({'data':res})
    return JsonResponse({})

# =============================2. Create New Organism based Taxonomoy==============================#
# OrgForm=modelform_factory(Organisms, exclude=["Strain_Type"])
def newOrgnisms(req):
   
    risk=Dictionaries.objects.filter(Dictionary_ID='Risk_Group')
    org_strain=Dictionaries.objects.filter(Dictionary_ID='unit_conversion')
    strains=str(org_strain[0]).replace("[","").replace("]","").replace('\'','')
    strain_test=strains.split(',')
    print(strain_test)
    a=strain_test[0]
    b=strain_test[1]
    
    if req.method=='POST':
        setTaxo=req.POST.get('setTaxo')
        print("setTaxo is type of ")
        print(setTaxo)
        newTaxo=get_object_or_404(Taxonomy, Organism_Name=setTaxo)
        # form=OrgForm(req.POST)
        strain=req.POST.get('strain')
        riskgroup=req.POST.get('risk')
        riskgroup=ast.literal_eval(riskgroup)
        # print(f'{riskgroup} type is: {type(riskgroup)}')
        try:
            riskg=get_object_or_404(Dictionaries, Dict_Value=riskgroup)
        except Exception as err:
            print(err)
            riskg=get_object_or_404(Dictionaries, pk=2)
        print(riskg.Dict_Value)

        try:
            newobj=Organisms.objects.create(Organism_Name=newTaxo, Strain_Type=strain, Risk_Group=riskg)
            newobj.save()
            print('saved!')
            return redirect("/home2")
            # if form.is_valid():
            #     instance=form.save()
            #     instance.save()
            #     obj=Organisms.objects.filter(Organism_Name='testSpeice')[0]
            #     obj.Strain_Type=strain
            #     obj.save(update_fields=['Strain_Type'])
            #     print(f'{obj.Strain_Type} saved')
                
            #     return redirect("org_create")
        except Exception as err:
            print(err)
    # else:
        # form=OrgForm()
    return render(req, 'aa_chem/orgCreate2.html', { 'a':a, 'b':b, 'risks':risk}) #'form':form,

# ==============================List View ===============================================================
# class OrgListView(ListView):
#     model=Organisms
#     paginate_by = 24
#     fields='__all__'
#     template_name = 'aa_chem/orgList.html'
    
#     def get_context_data(self, **kwargs):
#         for filename in os.listdir("static/images"):
#                 file_path=os.path.join("static/images", filename)
#                 try:
#                     os.unlink(file_path)
#                     print("removed!")
#                 except Exception as err:
#                     print(err)

#         context=super().get_context_data(**kwargs)
#         context["objects"]=self.model.objects.all()
#         objects_all=[object_ for object_ in context["objects"]]
#         p=Paginator(objects_all, 24)
        
#         for object_ in p.get_page(self.request.GET.get('page')):
       
#             # if exists(f"static/images/{object_.id}.png"):
#             #     print("file exists")
#             #     return context
#             # else:

#             if object_.id%2==0:
#                 m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
#             else:
#                 m=Chem.MolFromSmiles('CCCC[C@@H]1NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](N)Cc2ccccc2)C(C)C)CCC(=O)NCCCC[C@@H](C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc2c[nH]cn2)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCCN=C(N)N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H](C(=O)C(N)=O)[C@@H](C)CC)[C@@H](C)CC)NC(=O)[C@H](C)NC(=O)[C@H](CCCN=C(N)N)NC(=O)[C@H](C)NC1=O')
#             molecule_to_svg(m, object_.id)
            
#         return context

# class OrgTableView(OrgListView):
  
#     template_name = 'aa_chem/orgTable.html'

   



# # ==============================Function View===========================================================#


# # @staff_member_required
# # @user_passes_test(lambda u: u.is_staff)
# @permission_required('app.change_groupfilter')
# # @login_required
# def importCSV(req):  

#     if req.method=='POST':
#         try:
#             fdata=[req.POST.get('fdata'),]
#             management.call_command('processCsvmultiple', fpath=fdata)
#             print("done!")
#             print(req.body)
#         except Exception as err:
#             print(err)
   
#     return redirect(req.META['HTTP_REFERER'])

#======================================================Export Data Views Function==================================================#  
# # @user_passes_test(lambda u: u.is_admin)
# @permission_required('importdata')
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
   