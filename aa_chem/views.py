import os
from rdkit import Chem
from django_filters.views import FilterView

# from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib import messages
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction, IntegrityError
# from django.forms import modelform_factory
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView
from aa_chem.models import  Organisms, Taxonomy
from aa_chem.utils import  querysetToChoiseList_Dictionaries, searchbar_02, MySearchbar02, MySearchbar03
from app.models import Dictionaries
from .forms import CreateOrganism_form, UpdateOrganism_form, Taxonomy_form
# from coadd_web.settings import Strain_Type_choices




# # =======================================Taxonomy Read Create Update Delete View=============================================================================#

# Taxonomy Card View in Chem Homepage===============Read=================================================
def get_objects(req, model, model_field):
    req=req
    model=model
    model_field=model_field 
    objects=searchbar_02(req, model, model_field)
    p=Paginator(objects, 24)
    pag_num = req.GET.get('page')
    pag_obj=p.get_page(pag_num)
    
    context={
        'pag_obj':pag_obj,       
    }
    return context    



@login_required
def home(req): 
    context={}
    search_filter=MySearchbar02(req.GET, queryset=Taxonomy.objects.all())
    context['filter']=search_filter
    p=Paginator(search_filter.qs, 24)
    page_num=req.GET.get('page')
    page_obj=p.get_page(page_num)

    context['page_obj']=page_obj
    return render(req, 'aa_chem/chem.html', context)


# ==========List View================================Read===========================================
class TaxoListView(ListView):
    model=Taxonomy
    paginate_by = 24
    fields='__all__'
    template_name = 'aa_chem/readForm/Taxonomy_list.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.filter(astatus__gte=0)      
        return context

# ====================================================Create===========================================
def TaxoCreate(req):
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
def TaxoUpdate(req, pk):
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
# @user_passes_test(lambda u: u.is_superuser)
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

# ====================================================CREATE==========================================#
    # ==============Step1. Ajax Call search Taxonomy(for all models using Taxonomy as ForeignKey)=================#

    #  ===================refer to utils.py function searchbar_01==========================================

    # =============================step 2. Create new record by form===================#
def createOrgnisms(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''
    Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type']) # 
    kwargs={}
    kwargs['user']=req.user 
    if req.method=='POST':
        form=CreateOrganism_form(Strain_Type_choices,  req.POST,)
        Strain_Type_list=req.POST.getlist('Strain_Type')
        
        try:
            if form.is_valid():
                print("form is valid")  
                Organism_Name=req.POST.get('searchbar_01')
                print(f"request.Post.get {Organism_Name}")
                Organism = get_object_or_404(Taxonomy, Organism_Name=Organism_Name)
                form.get_object(Organism_Name) 
                instance=form.save(commit=False)
                print(f'instance.Organism_Name is {instance.Organism_Name} ')
                try:
                    with transaction.atomic():
                        instance.Organism_Name=Organism
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
        form=CreateOrganism_form(Strain_Type_choices,)
 
    return render(req, 'aa_chem/createForm/Organism.html', { 'form':form, 'Strain_Type':Strain_Type_choices})


#===============================================================Organism detail========================================================================================

def detailOrganism(req, pk):
    object_=get_object_or_404(Organisms, Organism_ID=pk)
    context={'Organism': object_}
    return render(req, "aa_chem/readForm/Organism_detail.html", context)

#======================================================================Update Organism=================================================================================

# @transaction.atomic(using='drugs_db')
def updateOrganism(req, pk):
    Strain_Type_choices=querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type'])
    object_=get_object_or_404(Organisms, Organism_ID=pk)
    original_Organism_Name=object_.Organism_Name
    

    kwargs={}
    kwargs['user']=req.user
    #This can be minimized when all organism have classes... 
    if object_.Organism_Name.Class:
        Organism_Class=object_.Organism_Name.Class.Dict_Value
        original_class=Organism_Class
    else:
        Organism_Class="No Class"
        original_class="No Class"

    if req.method=='POST':
        
        try:
            with transaction.atomic(using='drugs_db'):
                obj = Organisms.objects.select_for_update().get(Organism_ID=pk) 
                form=UpdateOrganism_form(Strain_Type_choices, req.POST, instance=obj)     
                if  req.POST.get('searchbar_01'):
                    Organism_Name=req.POST.get('searchbar_01')
                    print(f'http Request Organism_name is {Organism_Name}')

                    form.clean_organismName(Organism_Name, original_class)
                   
                else:
                    form.get_object(original_Organism_Name)
                        
                if form.is_valid():               
                    # If Update Organism_Name============================   
                    print("form is valid")             
                    instance=form.save(commit=False)
                    instance.Organism_Name=get_object_or_404(Taxonomy, Organism_Name=original_Organism_Name)  # here is a bug need to fix! 
              
                    instance.save(**kwargs)
                    print('save updated')
                   
                    return redirect(req.META['HTTP_REFERER'])

                else:
                    messages.warning(req, f"Form not Valid!{form.errors}")
                    return HttpResponse(status=204)
                   
        except Exception as err:
            messages.warning(req, f'{err} is the exception error')
            return redirect(req.META['HTTP_REFERER']) 
    

    else:
        form=UpdateOrganism_form(Strain_Type_choices, instance=object_)

    context={
        "form":form,
        "Organism":object_,
        "Class":Organism_Class
    }
   
    return render(req, "aa_chem/updateForm/Organism.html", context)

# ==============================Delete  ===============================================================
# @user_passes_test(lambda u: u.is_superuser) #login_url='/redirect/to/somewhere'
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
class OrganismListView(ListView):
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

class OrganismCardView(ListView):
    model=Organisms   
    template_name = 'aa_chem/readForm/Organism_card.html'
    paginate_by=3
    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context['filter']=MySearchbar03(self.request.GET, queryset=self.get_queryset())
        return context

    def get_queryset(self):
        queryset=super().get_queryset()
        return MySearchbar03(self.request.GET, queryset=queryset).qs

  
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
   