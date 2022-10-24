# import re
# import csv
import os
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
    else:
        result=Taxonomy.objects.all()
        
    objects_all=result
    p=Paginator(objects_all, 24)
    page_number = req.GET.get('page')
    page_obj=p.get_page(page_number)
    for object_ in page_obj:
        m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
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

            m=Chem.MolFromSmiles('Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1')
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

# =============================2. Create New Organism based Taxonomy==============================#

def newOrgnisms(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''
    strain_type=Dictionaries.objects.filter(Dictionary_ID='Strain_Type') #===multi choice

    # ===================Dictionary Foreign Key===================
    risk=Dictionaries.objects.filter(Dictionary_ID='Risk_Group')
    # pathogen=Dictionaries.objects.filter(Dictionary_ID='unit_conversion')
    # biolAppr =Dictionaries.objects.filter(Dictionary_ID='unit_conversion')
    # mat=Dictionaries.objects.filter(Dictionary_ID='unit_conversion')
    # oxyPref=Dictionaries.objects.filter(Dictionary_ID='unit_conversion')
    # ===================Dictionary Foreign Key===================

    #=====================normal fields==========================
    # Organism_Desc= models.CharField
    # Strain_ID= models.CharField
    # Strain_Code= models.CharField
    # Strain_Desc= models.CharField
    # Strain_Notes= models.CharField
    # Strain_Tissue= models.CharField    
    # Sequence = models.CharField
    # Sequence_Link = models.CharField
    # Tax_ID = models.IntegerField
    # Import_Permit = models.CharField
    # Special_Precaution = models.CharField
    # Lab_Restriction = models.CharField
    # MTA_Document = models.CharField
    # Atmosphere_Pref = models.CharField
    # Nutrient_Pref = models.CharFiel
    # Biofilm_Pref = models.CharField

    # Retreive Values for each column========================
    if req.method=='POST':
        setTaxo=req.POST.get('setTaxo')
        newTaxo=get_object_or_404(Taxonomy, Organism_Name=setTaxo)
        strain=req.POST.getlist('strain')
        print(type(strain))

        riskgroup=req.POST.get('risk')
        try:
            riskg=get_object_or_404(Dictionaries, Dict_Value=riskgroup)
            
        except Exception as err:
            print(err)
            riskg=None
        

        try:
            newobj=Organisms.objects.create(Organism_Name=newTaxo, Strain_Type=strain, Risk_Group=riskg)
            newobj.save()
            print('saved!')
            return redirect("/")
        
        except Exception as err:
            print(err)
 
    return render(req, 'aa_chem/orgCreate2.html', { 'strains':strain_type, 'risks':risk}) #'form':form,







# ==============================List View ===============================================================
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

# class OrgTableView(OrgListView):
  
#     template_name = 'aa_chem/orgTable.html'

   



# # ==============================Import Excel files===========================================================#
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage

def import_excel_taxo(req):
    print('importing....')
    try:
        if req.method=='POST' and req.FILES['myfile']:
            myfile=req.FILES['myfile']
            fs=FileSystemStorage()
            filename=fs.save(myfile.name, myfile)
            uploaded_file_url=fs.url(filename)
            excel_file=uploaded_file_url
            print(excel_file)
            exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
            print(type(exmpexceldata))
            dbframe=exmpexceldata
            for dbframe in dbframe.itertuples():
                class_fkey=Dictionaries.objects.filter(Dict_Desc=dbframe.ORGANISM_CLASS)
                print(class_fkey)
                division_fkey=Dictionaries.objects.filter(Dict_Desc=dbframe.DIVISION)
                linea=dbframe.LINEAGE.split(",")
                print(division_fkey)
                # fromdata_time_obj=dt.datetime.strptime(dbframe.DOB, '%d-%m-%Y')
                obj, created=Taxonomy.objects.get_or_create(Organism_Name=dbframe.ORGANISM_NAME, Other_Names=dbframe.ORGANISM_NAME_OTHER, Code=dbframe.ORGANISM_CODE, 
                    Class=class_fkey[0], Tax_ID=dbframe.TAX_ID,Parent_Tax_ID=dbframe.PARENT_TAX_ID, 
                    Tax_Rank=dbframe.TAX_RANK, Division=division_fkey[0], Lineage=linea
                    )
                print(type(obj))
                # obj.save()
            
            return render(req, 'aa_chem/importexcel_taxo.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/importexcel_taxo.html', {})

def import_excel_dict(req):
    print('importing....')
    try:
        if req.method=='POST' and req.FILES['myfile']:
            myfile=req.FILES['myfile']
            fs=FileSystemStorage()
            filename=fs.save(myfile.name, myfile)
            uploaded_file_url=fs.url(filename)
            excel_file=uploaded_file_url
            print(excel_file)
            exmpexceldata=pd.read_csv("."+excel_file, encoding='utf-8')
            print(type(exmpexceldata))
            dbframe=exmpexceldata
            for dbframe in dbframe.itertuples():
                print("iam here")             
                # int_order=int(0+dbframe.Dict_View_Order)
                # fromdata_time_obj=dt.datetime.strptime(dbframe.DOB, '%d-%m-%Y')
                obj, created=Dictionaries.objects.get_or_create(Dictionary_ID=dbframe.Dictionary_ID, Dictionary_Class=dbframe.Dictionary_Class, Dict_Value=dbframe.Dict_Value, Dict_Desc =dbframe.Dict_Desc, Dict_Value_Type =dbframe.Dict_Value_Type, Dict_View_Order =dbframe.Dict_View_Order
                    )
                print(type(obj))
                # obj.save()
            
            return render(req, 'aa_chem/importexcel_dict.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/importexcel_dict.html', {})

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
   