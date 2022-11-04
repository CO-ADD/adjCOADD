import os
from rdkit import Chem
from django_filters.views import FilterView

from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction, IntegrityError
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView

from aa_chem.models import  Organisms, Taxonomy
from aa_chem.utils import  querysetToChoiseList_Dictionaries, MySearchbar02, MySearchbar03
from app.models import Dictionaries
from .forms import CreateOrganism_form, UpdateOrganism_form, Taxonomy_form
from django.core.exceptions import ValidationError


# # =======================================Taxonomy Read Create Update Delete View=============================================================================#

# Taxonomy Card View in Chem Homepage===============Read=================================================
class TaxonomyCardView(LoginRequiredMixin, ListView):
    model=Taxonomy  
    template_name = 'aa_chem/readForm/Taxonomy_card.html' 
    paginate_by=24

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context['filter']=MySearchbar02(self.request.GET, queryset=self.get_queryset())
        return context

    def get_queryset(self):
        qs=super().get_queryset()
        return MySearchbar02(self.request.GET, queryset=qs).qs

# ==========List View================================Read===========================================
class TaxonomyListView(TaxonomyCardView):
    template_name = 'aa_chem/readForm/Taxonomy_list.html'

# ====================================================Create===========================================
@login_required
@user_passes_test(lambda u: u.is_staff) 
def createTaxonomy(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Taxonomy_form
    if req.method=='POST':
        form=Taxonomy_form(req.POST)
        if form.is_valid():
            print("form is valid")   
            instance=form.save(commit=False)
            instance.save(**kwargs)
            print("saved")
            return redirect("/")
        else:
            print(form.errors)
    return render(req, 'aa_chem/createForm/Taxonomy.html', {'form':form})
    
# ====================================================Update===========================================
@login_required
@user_passes_test(lambda u: u.is_staff) 
def updateTaxonomy(req, pk):
    object_=get_object_or_404(Taxonomy, Organism_Name=pk)
    kwargs={}
    kwargs['user']=req.user 
    form=Taxonomy_form(instance=object_)
    if req.method=='POST':
        form=Taxonomy_form(req.POST, instance=object_)
        if form.is_valid():
            print("form is valid")   
            instance=form.save(commit=False)        
            instance.save(**kwargs)
            print("saved")
            return redirect("/")
        else:
            print(form.errors)
    return render(req, 'aa_chem/updateForm/Taxonomy.html', {'form':form, 'object':object_})

# ====================================================Delete===========================================
@user_passes_test(lambda u: u.is_superuser)
def deleteTaxonomy(req, pk):
    kwargs={}
    kwargs['user']=req.user 
    print('deleting view')
    object_=get_object_or_404(Taxonomy, Organism_Name=pk)
    try:
        print(object_.Organism_Name)
        object_.delete(**kwargs)
        print("deleted")
    except Exception as err:
        print(err)
    return redirect("/")

# # ========================================Organisms CREATE READ UPDATE DELETE View==============================================#

# ======================================================================CREATE==========================================#
    # ==============Step1. Ajax Call search Taxonomy(for all models using Taxonomy as ForeignKey)=================#

    #  ===================refer to utils.py function searchbar_01==========================================

    # =============================step 2. Create new record by form===================#
@login_required
@user_passes_test(lambda u: u.is_staff) 
def createOrgnisms(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''
    
    kwargs={}
    kwargs['user']=req.user 
    if req.method=='POST':
        Organism_Name=req.POST.get('searchbar_01')
        Strain_Type_list=req.POST.getlist('Strain_Type')
        form=CreateOrganism_form(Organism_Name, req.POST)
        print(f"request.Post.get {Organism_Name}")     
        
        try:
            if form.is_valid():
                print("form is valid")  
                try:
                    with transaction.atomic():
                        instance=form.save(commit=False)                  
                        instance.save(**kwargs)
                        print("saved")
                        return redirect("org_list")

                except IntegrityError as err:
                        messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                        return redirect(req.META['HTTP_REFERER'])                
            else:
                print(f'something wrong...{form.errors}')
                return redirect(req.META['HTTP_REFERER'])      
        except Exception as err:
            print(f'error is {form.errors} with {err}')
            return redirect(req.META['HTTP_REFERER'])

    else:
        form=CreateOrganism_form()
    return render(req, 'aa_chem/createForm/Organism.html', { 'form':form, }) 


#=====================================================================Organism detail========================================================================================
@login_required
def detailOrganism(req, pk):
    object_=get_object_or_404(Organisms, Organism_ID=pk)
    context={'Organism': object_}
    return render(req, "aa_chem/readForm/Organism_detail.html", context)

#======================================================================Update Organism=================================================================================
@login_required
@user_passes_test(lambda u: u.is_staff) 
def updateOrganism(req, pk):
    object_=get_object_or_404(Organisms, Organism_ID=pk)
    kwargs={}
    kwargs['user']=req.user
    #This can be minimized when all organism have classes... ----------------
    if object_.Organism_Name.Class:
        Organism_Class_str=object_.Organism_Name.Class.Dict_Value
    else:
        Organism_Class_str="No Class"
    #-------------------------------------------------------------------------
    if req.method=='POST':

        try:
            with transaction.atomic(using='drugs_db'):        # testing!
                obj = Organisms.objects.select_for_update().get(Organism_ID=pk)
                #------------------------If update Organism Name-----------------------------------
                if  req.POST.get('searchbar_01'):
                    Organism_Name_str=req.POST.get('searchbar_01')
                    Organism_new_obj=get_object_or_404(Taxonomy, Organism_Name=Organism_Name_str)
                    form=UpdateOrganism_form(Organism_Name, req.POST, instance=obj)
                    print('form created')
                    #-----------------------Not allow to update name in different class--------
                    if Organism_new_obj.Class.Dict_Value and Organism_new_obj.Class.Dict_Value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                    #-----------------------Not allow to update name in different class--------
                #------------------------If update Organism Name-----------------------------------
                else:
                    form=UpdateOrganism_form(object_.Organism_Name.Organism_Name, req.POST, instance=obj) 
                try:
                    if form.is_valid():                  
                        instance=form.save(commit=False)
                        instance.save(**kwargs)
                        print('save updated')                 
                        return redirect(req.META['HTTP_REFERER'])
                except Exception as err:
                    print(err)
                   
        except Exception as err:
            messages.warning(req, f'Update failed due to {err} error')
            return redirect(req.META['HTTP_REFERER'])
  
    else:
        form=UpdateOrganism_form(instance=object_)

    context={
        "form":form,
        "Organism":object_,
        "Class":Organism_Class_str
    }
   
    return render(req, "aa_chem/updateForm/Organism.html", context)

# ==============================Delete  ===============================================================
@user_passes_test(lambda u: u.is_superuser) #login_url='/redirect/to/somewhere'
def deleteOrganism(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(Organisms, Organism_ID=pk)
    try:      
        object_.delete(**kwargs)
        print("deleted")
    except Exception as err:
        print(err)
    return redirect("/")
    

# ==============================List View ===============================================================
class OrganismListView(LoginRequiredMixin, ListView):
    model=Organisms  
    template_name = 'aa_chem/readForm/Organism_list.html' 
    paginate_by=3

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context['filter']=MySearchbar03(self.request.GET, queryset=self.get_queryset())
        return context

    def get_queryset(self):
        queryset=super().get_queryset()
        return MySearchbar03(self.request.GET, queryset=queryset).qs

class OrganismCardView(OrganismListView):

    template_name = 'aa_chem/readForm/Organism_card.html'
  
    

  
# # ==============================Import Excel files===========================================================#
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage
@login_required
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
                class_fkey=Dictionaries.objects.filter(Dict_Value=dbframe.ORGANISM_CLASS)
                if class_fkey:
                    class_fkey=class_fkey[0]
                else:
                    class_fkey=None
                print(class_fkey)
                division_fkey=Dictionaries.objects.filter(Dict_Value=dbframe.DIVISION)
                if division_fkey:
                    division_fkey=division_fkey[0]
                else:
                    division_fkey=None
                linea=str(dbframe.LINEAGE).split(";")
                print(division_fkey)
                # fromdata_time_obj=dt.datetime.strptime(dbframe.DOB, '%d-%m-%Y')
                try:
                    obj, created=Taxonomy.objects.get_or_create(Organism_Name=dbframe.ORGANISM_NAME, Other_Names=dbframe.ORGANISM_NAME_OTHER, Code=dbframe.ORGANISM_CODE, 
                        Class=class_fkey, Tax_ID=dbframe.TAX_ID, Parent_Tax_ID=dbframe.PARENT_TAX_ID, 
                        Tax_Rank=dbframe.TAX_RANK, Division=division_fkey, Lineage=linea
                        )
                except Exception as err:
                    print(err)
                # obj.save()
            
            return render(req, 'aa_chem/createForm/importDataForm/importexcel_taxo.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/createForm/importDataForm/importexcel_taxo.html', {})
#=======================================================================================================
@login_required
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
                obj, created=Dictionaries.objects.get_or_create(Dictionary_Class=dbframe.Class, Dict_Value=dbframe.Term, Dict_Desc =dbframe.Name)
                print(type(obj))
                # obj.save()
            
            return render(req, 'aa_chem/createForm/importDataForm/importexcel_dict.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'aa_chem/createForm/importDataForm/importexcel_dict.html', {})

#======================================================Export Data Views Function==================================================#  
# # @user_passes_test(lambda u: u.is_admin)
# @permission_required('importdata')
@login_required
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
   