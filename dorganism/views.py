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

from .models import  Organism, Taxonomy, Organism_Batch, OrgBatch_Stock
from .utils import  Organismfilter, Taxonomyfilter, Batchfilter
from apputil.models import Dictionary, ApplicationUser
from apputil.views import permission_not_granted
from .forms import CreateOrganism_form, UpdateOrganism_form, Taxonomy_form, Batch_form, Batchupdate_form, Stock_form

#  #####################Django Filter View#################
# Base Class for all models list/card view
class FilteredListView(ListView):
    filterset_class = None
    paginate_by=50

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
        context['fields']=self.model.get_fields()
        print(context['fields'])
        return context

    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)
        return paginate_by


# #############################TAXONOMY View############################################
# ==========List View================================Read===========================================
class TaxonomyListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Taxonomy  
    template_name = 'dorganism/readForm/Taxonomy_list.html' 
    filterset_class=Taxonomyfilter

 
class TaxonomyCardView(TaxonomyListView):
    template_name = 'dorganism/readForm/Taxonomy_card.html'
# ===========Detail View=============================Read============================================
@login_required
def detailTaxonomy(req, slug=None):
    context={}
    object_=get_object_or_404(Taxonomy, urlname=slug)
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
def updateTaxonomy(req, slug=None):
    object_=get_object_or_404(Taxonomy, urlname=slug)
    kwargs={}
    kwargs['user']=req.user 
    form=Taxonomy_form(instance=object_)
    if req.method=='POST':
        form=Taxonomy_form(req.POST, instance=object_)
        if form.is_valid():
            instance=form.save(commit=False)        
            instance.save(**kwargs)
            print("saved")
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'dorganism/updateForm/Taxonomy_u.html', {'form':form, 'object':object_})

# ====================================================Delete===========================================
@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted')
def deleteTaxonomy(req, slug=None):
    kwargs={}
    kwargs['user']=req.user 
    object_=get_object_or_404(Taxonomy, urlname=slug)
    try:
        if req.method=='POST':
            object_.delete(**kwargs)
    except Exception as err:
        print(err) 
    return redirect("taxo_card")

############################################### ORGANISM View ###########################################
# ==============================List View ================================
class OrganismListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism  
    template_name = 'dorganism/readForm/Organism_list.html'
    filterset_class=Organismfilter
    
class OrganismCardView(OrganismListView):
    template_name = 'dorganism/readForm/Organism_card.html'

# ======================================================================CREATE==========================================#
    # ==Step1. Ajax Call(def search_organism in utils) search Taxonomy(for all models using Taxonomy as ForeignKey)=====#
    # =============================step 2. Create new record by form===================#
# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createOrganism(req):
    '''
    Function View Create new Organism table row with foreignkey: Taxonomy and Dictionary. 
    '''  
    kwargs={}
    kwargs['user']=req.user

    form=CreateOrganism_form(req.user)
    if req.method=='POST':
        Organism_Name=req.POST.get('search_organism')
        form=CreateOrganism_form( req.user, Organism_Name, req.POST,)
        print(f"request.Post.get {Organism_Name}")     
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
    context["batch_obj"]=Organism_Batch.objects.filter(organism_id=object_.organism_id, astatus__gte=0)
    context["batch_fields"]=Organism_Batch.get_fields()



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
   
# #############################BATCH View############################################
# ==========List View================================Read===========================================
class BatchCardView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Organism_Batch 
    template_name = 'dorganism/readForm/OrganismBatch_card.html' 
    filterset_class=Batchfilter

    
# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createBatch(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Batch_form(req.user)

    if req.method=='POST':
        Organism_Id=req.POST.get('search_organism')
        form=Batch_form(req.user, Organism_Id, req.POST)
        if form.is_valid():
            print("form is valid")  
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    # print(instance.organism_id)                 
                    instance.save(**kwargs)
                    print("new Batch saved--view info")
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
        

    return render(req, 'dorganism/createForm/Batch_c.html', { 'form':form, }) 

from django.http import QueryDict

@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateBatch(req, pk):
    print(req.method)
    object_=get_object_or_404(Organism_Batch, orgbatch_id=pk)
    kwargs={}
    kwargs['user']=req.user
   
    form=Batchupdate_form(req.user, instance=object_)
   
    context={
        "form":form,
        "object":object_,
    }
    #-------------------------------------------------------------------------
    
    if req.method=='PUT':
        qd=QueryDict(req.body).dict()
        print(qd)
        object_batch=get_object_or_404(Organism_Batch, orgbatch_id=qd["orgbatch_id"])
        form=Batchupdate_form(req.user, data=qd, instance=object_batch, )
        print(qd)
        
        if form.is_valid():
            kwargs={}
            kwargs['user']=req.user                  
            instance=form.save(commit=False)
            instance.save(**kwargs)
            context={
                "object_batch":object_batch,
                'object':object_batch  # this object refer to the same entry of object_batch
            }
            return render(req, "dorganism/readForm/Batch_tr.html", context)
            # return render(req, "dorganism/updateForm/Batch_u.html", context)  
                
   
    return render(req, "dorganism/updateForm/Batch_u.html", context)

@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted') 
def deleteBatch(req, pk):
    kwargs={}
    kwargs['user']=req.user
    print(f'batchID {pk}')
    object_=get_object_or_404(Organism_Batch, orgbatch_id=pk)
    try:
        object_.delete(**kwargs)
        print("deleted")
            
    except Exception as err:
        print(err)
    return redirect('/')


############################################### Stock View ###########################################
# class StockListView(LoginRequiredMixin, ListView):
#     login_url = '/'
#     model=OrgBatch_Stock
#     template_name = 'dorganism/readForm/OrganismStock_list.html' 

#     def get_context_data(self, **kwargs):
#         context = super().get_context_data(**kwargs)
#         project = OrgBatch_Stock.objects.filter(orgbatch_id=self.kwargs.get('pk'), astatus__gte=0)
        
#         context["object_list"]=project
#         context["stock_fields"]=OrgBatch_Stock.get_fields()
#         return context
@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted') 
def stockList(req, pk):
    # if req.headers.get('x-requested-with') == 'XMLHttpRequest':
    res=None
    if req.method == 'GET':
        batch_id=req.GET.get('Batch_id')
        print(f"StockList with ID = {batch_id}")
        object_=Organism_Batch.objects.get(orgbatch_id=batch_id)
        qs=OrgBatch_Stock.objects.filter(orgbatch_id=object_, astatus__gte=0)
        if len(qs)>0:
            data=[]
            for i in qs:
                item={
                    "stock_id":i.pk,
                    "stock_created":i.n_created,
                    "stock_left":i.n_left,
                    "stock_note":i.stock_note,
                    "stock_type":i.stock_type.dict_value,
                    "stock_date":i.stock_date,
                }
                data.append(item)
            res=data
            print(res)
        else:
            res='No Data'
        return JsonResponse({'data':res})
    return JsonResponse({})

    

# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createStock(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Stock_form()

    if req.method=='POST':
        
        form=Stock_form(req.POST)
        if form.is_valid():
            print("form is valid")  
            try:
                with transaction.atomic(using='dorganism'):
                    instance=form.save(commit=False) 
                    print("form save")                 
                    instance.save(**kwargs)
                    print("saved--view info")
                    return redirect(req.META['HTTP_REFERER']) 

            except IntegrityError as err:
                    messages.error(req, f'IntegrityError {err} happens, record may be existed!')
                    return redirect(req.META['HTTP_REFERER'])                
        else:
            print(f'something wrong...{form.errors}')
            return redirect(req.META['HTTP_REFERER'])      
        

    return render(req, 'dorganism/createForm/Stock_c.html', { 'form':form, }) 


@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateStock(req, pk):
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    print(object_)
    kwargs={}
    kwargs['user']=req.user
   
    form=Stock_form(instance=object_)
    #-------------------------------------------------------------------------
    if req.method=='POST':
        form=Stock_form(req.POST, instance=object_)
        if "cancel" in req.POST:
            return redirect(req.META['HTTP_REFERER'])
        else:

            try:
                with transaction.atomic(using='dorganism'):        # testing!
                    obj = OrgBatch_Stock.objects.select_for_update().get(pk=pk)
                    try:
                        if form.is_valid():                  
                            instance=form.save(commit=False)
                            instance.save(**kwargs)
                            print('updated')
                            return redirect(req.META['HTTP_REFERER'])
                    except Exception as err:
                        print(f'form erroro is {form.errors} and error {err}')
                   
            except Exception as err:
                messages.warning(req, f'Update failed due to {err} error')
                return redirect(req.META['HTTP_REFERER'])
    context={
        "form":form,
        "object":object_,
    }
   
    return render(req, "dorganism/updateForm/Stock_u.html", context)

@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted') 
def deleteStock(req, pk):
    kwargs={}
    kwargs['user']=req.user
    object_=get_object_or_404(OrgBatch_Stock, pk=pk)
    context={'object':object_}
    # if request.method == 'GET':
    if req.method=='POST':
        object_.delete(**kwargs)
        print("deleted")
        return redirect(req.META['HTTP_REFERER'])
    return render(req, "dorganism/deleteForm/Stock_del.html", context)







############################################### Import CSV View ###########################################
import pandas as pd
from django.conf import settings
from django.core.files.storage import FileSystemStorage
@login_required
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted')
@transaction.atomic
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
            
            return render(req, 'dorganism/createForm/importDataForm/importexcel.html', {'uploaded_file_url': uploaded_file_url})
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


############################################### Export CSV View ###########################################
import csv
import datetime 

@login_required
@user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
def exportCSV(req):
    queryset=Taxonomy.objects.all()
    query= Taxonomyfilter(req.GET, queryset=queryset).qs
    response = HttpResponse(content_type='text/csv')
    file_name = "fltred_loaction_data" + str(datetime.date.today()) + ".csv"

    writer = csv.writer(response)
    writer.writerow(['name', 'class','lineage',])
    for i in query.values_list('organism_name','org_class', 'lineage',):
        writer.writerow(i)
    response['Content-Disposition'] = 'attachment; filename = "' + file_name + '"'
    return response