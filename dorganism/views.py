import os
from rdkit import Chem
from django_filters.views import FilterView

from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import transaction, IntegrityError
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.urls import reverse_lazy
from django.views.decorators.csrf import csrf_exempt, csrf_protect
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView, TemplateView
from django.utils.functional import SimpleLazyObject

from .models import  Organism, Taxonomy
from .utils import  Organismfilter, Taxonomyfilter
from apputil.models import Dictionary, ApplicationUser
from apputil.views import permission_not_granted
from .forms import CreateOrganism_form, UpdateOrganism_form, Taxonomy_form

#  #####################Django Filter View#################
class FilteredListView(ListView):
    filterset_class = None
    paginate_by=10

    def get_queryset(self):
        # Get the queryset however you usually would.  For example:
        queryset = super().get_queryset()
        # Then use the query parameters and the queryset to
        # instantiate a filterset and save it as an attribute
        # on the view instance for later.
        self.filterset = self.filterset_class(self.request.GET, queryset=queryset)
        # Return the filtered queryset
        return self.filterset.qs.distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Pass the filterset to the template - it provides the form.
        context['filter'] = self.filterset
        context['paginate_by']=self.get_paginate_by(self, **kwargs)
        return context

    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)
        return paginate_by


# #############################TAXONOMY View############################################
# ==========List View================================Read===========================================
class TaxonomyListView(LoginRequiredMixin, FilteredListView):
    model=Taxonomy  
    template_name = 'dorganism/readForm/Taxonomy_list.html' 
    filterset_class=Taxonomyfilter

 
class TaxonomyCardView(TaxonomyListView):
    template_name = 'dorganism/readForm/Taxonomy_card.html'
    paginate_by=24
 
 
# ===========Detail View=============================Read============================================
@login_required
def detailTaxonomy(req, pk):
    context={}
    object_=get_object_or_404(Taxonomy, organism_name=pk)
    context["object"]=object_
    context['form']=Taxonomy_form(instance=object_)
    return render(req, "dorganism/readForm/Taxonomy_detail.html", context)

# ====================================================Create===========================================
# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
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
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dorganism/createForm/Taxonomy_c.html', {'form':form})
    
# ====================================================Update in Form===========================================
@login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateTaxonomy(req, pk):
    object_=get_object_or_404(Taxonomy, organism_name=pk)
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
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'dorganism/updateForm/Taxonomy_u.html', {'form':form, 'object':object_})

# ====================================================Delete===========================================
@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted')
def deleteTaxonomy(req, pk):
    kwargs={}
    kwargs['user']=req.user 
    print('deleting view')
    object_=get_object_or_404(Taxonomy, organism_name=pk)
    try:
        if req.method=='POST':
            object_.delete(**kwargs)
            print("deleted")
    except Exception as err:
        print(err) 
    return redirect("taxo_card")

############################################### ORGANISM View ###########################################
# ==============================List View ================================
class OrganismListView(LoginRequiredMixin, FilteredListView):
    model=Organism  
    template_name = 'dorganism/readForm/Organism_list.html'
    filterset_class=Organismfilter
    paginate_by=10

    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)
        return paginate_by
    
class OrganismCardView(OrganismListView):
    template_name = 'dorganism/readForm/Organism_card.html'

# ======================================================================CREATE==========================================#
    # ==Step1. Ajax Call(def search_organism in utils) search Taxonomy(for all models using Taxonomy as ForeignKey)=====#
    # =============================step 2. Create new record by form===================#
@login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createOrganism(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user
    print(f"in view: {req.user}")

    form=CreateOrganism_form(req.user)
    if req.method=='POST':
        Organism_Name=req.POST.get('search_organism')
        Strain_Type_list=req.POST.getlist('strain_type')
        form=CreateOrganism_form( req.user, Organism_Name, req.POST,)
        print(f"request.Post.get {Organism_Name}")     
        
        try:
            if form.is_valid():
                print("form is valid")  
                try:
                    with transaction.atomic(using='dorganism'):
                        instance=form.save(commit=False) 
                        print("form save")                 
                        instance.save(**kwargs)
                        print("saved--view info")
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

    return render(req, 'dorganism/createForm/Organism_c.html', { 'form':form, }) 


#=========================================Organism detail table with updating in detail table========================================================================================
@login_required
def detailOrganism(req, pk):
    context={}
    object_=get_object_or_404(Organism, organism_id=pk)
    user=req.user
    form=UpdateOrganism_form(user,instance=object_)
   
    context["object"]=object_
    context["form"]=form

    return render(req, "dorganism/readForm/Organism_detail.html", context)

#======================================================================Update Organism=================================================================================
@login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateOrganism(req, pk):
    object_=get_object_or_404(Organism, organism_id=pk)
    kwargs={}
    kwargs['user']=req.user
    #This can be minimized when all organism have classes... ----------------
    if object_.organism_name.org_class:
        Organism_Class_str=object_.organism_name.org_class.dict_value
    else:
        Organism_Class_str="No Class"
    #-------------------------------------------------------------------------
    if req.method=='POST':

        try:
            with transaction.atomic(using='dorganism'):        # testing!
                obj = Organism.objects.select_for_update().get(organism_id=pk)
                #------------------------If update Organism Name-----------------------------------
                if  req.POST.get('search_organism'):
                    print(req.POST.get('search_organism'))
                    Organism_Name_str=req.POST.get('search_organism')
                    Organism_new_obj=get_object_or_404(Taxonomy, organism_name=Organism_Name_str)
                    form=UpdateOrganism_form(req.user, Organism_Name_str, req.POST, instance=obj)
                    print('form created')
                    #-----------------------Not allow to update name in different class--------
                    if Organism_new_obj.org_class.dict_value and Organism_new_obj.org_class.dict_value != Organism_Class_str:
                        raise ValidationError('Not the same Class')
                    #-----------------------Not allow to update name in different class--------
                #------------------------If update Organism Name-----------------------------------
                else:
                    form=UpdateOrganism_form(req.user, object_.organism_name, req.POST, instance=obj) 
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
        form=UpdateOrganism_form(req.user,instance=object_)

    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dorganism/updateForm/Organism_u.html", context)

# ==============================Delete  ===============================================================
@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted') 
def deleteOrganism(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(Organism, organism_id=pk)
    try:
        object_.delete(**kwargs)
        print("deleted")
            
    except Exception as err:
        print(err)
    return redirect('/')
   

# # ==============================Import Excel files===========================================================#
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage
@login_required
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
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
                class_fkey=Dictionary.objects.filter(dict_value=dbframe.ORGANISM_CLASS)
                if class_fkey:
                    class_fkey=class_fkey[0]
                else:
                    class_fkey=None
                print(class_fkey)
                division_fkey=Dictionary.objects.filter(dict_value=dbframe.DIVISION)
                if division_fkey:
                    division_fkey=division_fkey[0]
                else:
                    division_fkey=None
                linea=str(dbframe.LINEAGE).split(";")
                print(division_fkey)
                try:
                    obj, created=Taxonomy.objects.get_or_create(organism_name=dbframe.ORGANISM_NAME, other_names=dbframe.ORGANISM_NAME_OTHER, code=dbframe.ORGANISM_CODE, 
                        org_class=class_fkey, tax_id=dbframe.TAX_ID, parent_tax_id=dbframe.PARENT_TAX_ID, 
                        tax_rank=dbframe.TAX_RANK, division=division_fkey, lineage=linea, 
                        )
                except Exception as err:
                    print(err)
                # obj.save()
            
            return render(req, 'dorganism/createForm/importDataForm/importexcel_taxo.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {})
#=======================================================================================================
@login_required
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
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
                obj, created=Dictionary.objects.get_or_create(dict_class=dbframe.Class, dict_value=dbframe.Term, dict_desc =dbframe.Name, )
                print(type(obj))
          
            return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(f'import failed because {err}')
    return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {})



#==================================================================import Organism================================================
@login_required
@user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
def import_excel_organism(req):
    print('importing....')
    try:
        if req.method=='POST' and req.FILES['myfile']:
            myfile=req.FILES['myfile']
            fs=FileSystemStorage()
            filename=fs.save(myfile.name, myfile)
            uploaded_file_url=fs.url(filename)
            excel_file=uploaded_file_url
            print(excel_file)
            exmpexceldata=pd.read_csv("."+excel_file, )
            print(exmpexceldata.itertuples)
            dbframe=exmpexceldata
            for dbframe in dbframe.itertuples():
                taxID=int('0'+dbframe[22])
                screen_panel=dbframe[26].split(';')
                organism_fkey=Taxonomy.objects.filter(organism_name=dbframe[1])
                print(organism_fkey[0])   
                try:
                    obj, created=Organism.objects.get_or_create(organism_id=dbframe[0], organism_name=organism_fkey[0],  strain_id=dbframe[3], 
                                    strain_code=dbframe[5], strain_notes=dbframe[7], 
                                    strain_tissue=dbframe[25], strain_type=dbframe[4], sequence=dbframe[28], sequence_link=dbframe[29], 
                                    strain_panel=screen_panel, 
                                    tax_id =taxID,risk_group=dbframe[9], pathogen_group =dbframe[10],import_permit =dbframe[12],bio_approval =dbframe[23],special_precaution =dbframe[24],lab_restriction =dbframe[27],mta_document =dbframe[31],
                                    mta_status =dbframe[32],oxygen_pref =dbframe[13],atmosphere_pref ='containSpecialCHA', nutrient_pref =dbframe[15],biofilm_pref =dbframe[16], acreated_by=req.user )
                except Exception as err:
                    print(err)
                # obj.save()
            
            return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {'uploaded_file_url': uploaded_file_url})
    except Exception as err:
        print(err)
    return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {})


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
   