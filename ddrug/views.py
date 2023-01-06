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
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from django.views.generic import ListView, TemplateView
from django.utils.functional import SimpleLazyObject

from apputil.models import Dictionary, ApplicationUser
from apputil.utils import FilteredListView
from apputil.views import permission_not_granted
from adjcoadd.constants import *
from .models import  Drug
# from .utils import 
from .forms import Drug_form
   
          
# #############################Drug View############################################
# ==========List View================================Read===========================================
class DrugListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Drug  
    template_name = 'ddrug/drug/drug_list.html' 
    filterset_class=None
    model_fields=None

 
class DrugCardView(DrugListView):
    template_name = 'ddrug/drug/drug_card.html'

    
# ===========Detail View=============================Read============================================
@login_required
def detailDrug(req, pk):
    context={}
    object_=get_object_or_404(Drug,pk=pk)
    context["object"]=object_
    context['form']=Drug_form(instance=object_)
    return render(req, "ddrug/drug/drug_detail.html", context)

# ====================================================Create===========================================
# @login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def createDrug(req):
    kwargs={}
    kwargs['user']=req.user 
    form=Drug_form
    if req.method=='POST':
        form=Drug_form(req.POST)
        if form.is_valid():
            instance=form.save(commit=False)
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'ddrug/drug/drug_c.html', {'form':form})
    
# ====================================================Update in Form===========================================
@login_required
@user_passes_test(lambda u: u.has_permission('Write'), login_url='permission_not_granted') 
def updateDrug(req, pk):
    object_=get_object_or_404(Drug, pk=pk)
    kwargs={}
    kwargs['user']=req.user 
    form=Drug_form(instance=object_)
    if req.method=='POST':
        form=Drug_form(req.POST, instance=object_)
        if form.is_valid():
            instance=form.save(commit=False)        
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'ddrug/drug/drug_u.html', {'form':form, 'object':object_})

# ====================================================Delete===========================================
@user_passes_test(lambda u: u.has_permission('Delete'), login_url='permission_not_granted')
def deleteDrug(req, pk=pk):
    kwargs={}
    kwargs['user']=req.user 
    object_=get_object_or_404(Drug, pk=pk)
    try:
        if req.method=='POST':
            object_.delete(**kwargs)
    except Exception as err:
        print(err) 
    return redirect("drug_card")

# Create your views here.
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
#         'chose':Taxonomy.Choice_Dictionary   
       
#     }
  
#     return render(req, 'organism/chem.html', context)
