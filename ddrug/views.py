import json
import os
from rdkit import Chem
from django_rdkit.models import *
from django_rdkit.config import config
import pandas as pd
import numpy as np
from django.core.serializers.json import DjangoJSONEncoder
import logging
logger = logging.getLogger("django")

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render
from django.conf import settings

from adjcoadd.constants import *
from apputil.utils.filters_base import FilteredListView
from apputil.utils.api_class import API_ListView
from apputil.utils.views_base import SimplecreateView, SimpleupdateView
from adjcoadd.constants import *
from ddrug.models import  Drug, VITEK_AST, VITEK_Card, VITEK_ID, MIC_COADD, MIC_Pub, Breakpoint
from ddrug.utils.molecules import molecule_to_svg, get_mfp2_neighbors
from ddrug.forms import Drug_form, Drug_filter, VitekCard_Filter, VitekAST_Filter, VitekID_Filter,MIC_COADDfilter, MIC_Pubfilter, Breakpointfilter
from ddrug.serializers import Drug_Serializer, VITEK_ASTSerializer

# ===================================================================
@login_required   
def smartsQuery(req, pk):
    '''
    MAKE SUBSTRUCTURE QUERY 
    '''
    context={}
    object_=get_object_or_404(Drug, drug_id=pk)
    context["object"]=object_
    # get mol block for an object
    try:
        context["object_mol"]=Chem.MolToMolBlock(object_.smol)
    # convert object to JMSE regonized form
        m="\\n".join(context["object_mol"].split("\n")) 
        context["object_mol"]=m
    except Exception as err:
        logger.error(err)
        messages.error(req, f'{object_.pk} mol not exists or {err}')
        context["object_mol"]=''

    return render(req, "ddrug/drug/drug_detail_structure.html", context)

@login_required   
def iframe_url(req):
    context={}
    return render(req, "utils/ketcher/index.html")

@login_required   
def ketcher_test(req):
    context={}
    n=Chem.MolFromSmiles("CC(C)([C@@H]1C(O)=O)S[C@H]([C@@H]2NC([C@@H](c(cc3)ccc3O)N)=O)N1C2=O")
    m=Chem.MolToMolBlock(n)
    context["object_mol"]="\\n".join(m.split("\n"))
    return render(req, "utils/ketcher_test.html", context)

#==  Drug View =============================================================
#--  DrugList --------------------------------------------------------------
class DrugListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Drug  
    template_name = 'ddrug/drug/drug_list.html' 
    filterset_class=Drug_filter
    model_fields=model.HEADER_FIELDS

#--  DrugCard --------------------------------------------------------------
class DrugCardView(DrugListView):
    template_name = 'ddrug/drug/drug_card.html'
    def get_context_data(self, **kwargs):
        try:
            context = super().get_context_data(**kwargs)          
        # clearIMGfolder()
            for object_ in context["object_list"]:
                filepath=os.path.join(settings.STRUCTURE_FILES_DIR, f"{object_.pk}.svg") 
                if os.path.exists(filepath):
                    continue
                else:
                    m=object_.smol
                    try:
                        molecule_to_svg(m, object_.pk)
                    except Exception as err:
                        pass
                        # messages.error(self.request, f'**{object_.pk} mol may not exists**')
        except Exception as err:
            context={}
            messages.error(self.request, err)
        return context    

##
@login_required
def detailDrug(req, pk):
    context={}
    object_=get_object_or_404(Drug, pk=pk)
    smol_initial = Chem.MolToMolBlock(object_.smol) if object_.smol else None
    form=Drug_form(instance=object_, initial={"smol":smol_initial},)
    context["object"]=object_
    context["form"]=form
    context["Links"]=LinkList
    try:
        context["object_mol"]=Chem.MolToMolBlock(object_.smol)
        m="\\n".join(context["object_mol"].split("\n"))
        context["object_mol"]=m
    except Exception as err:
        context["object_mol"]=""
    return render(req, "ddrug/drug/drug_detail.html", context)

##
#--  DrugCreate --------------------------------------------------------------
class DrugCreateView(SimplecreateView):
    form_class=Drug_form
    template_name='ddrug/drug/drug_c.html'
    
##
#--  DrugUpdate --------------------------------------------------------------
class DrugUpdateView(SimpleupdateView):
    form_class=Drug_form
    template_name='ddrug/drug/drug_u.html'
    model=Drug


#=================================================================================================
# Vitek Data
#=================================================================================================

# -----------------------------------------------------------------
# VitekCard
# -----------------------------------------------------------------
class VitekCard_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_Card  
    template_name = 'ddrug/vitek_card/vitekcard_list.html' 
    filterset_class=VitekCard_Filter
    model_fields=model.HEADER_FIELDS
    #context_list=''
    
    # def get_context_data(self,  **kwargs):

    #     context =super().get_context_data( **kwargs)
    #     return context

# -----------------------------------------------------------------
# Vitek AST
# -----------------------------------------------------------------
class VitekAST_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_AST  
    template_name = 'ddrug/vitek_ast/vitekast_list.html' 
    filterset_class=VitekAST_Filter
    model_fields=model.HEADER_FIELDS
    #context_list=''

      
# -----------------------------------------------------------------
# Vitek ID
# -----------------------------------------------------------------
class VitekID_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=VITEK_ID 
    template_name = 'ddrug/vitek_id/vitekid_list.html' 
    filterset_class=VitekID_Filter
    model_fields=model.HEADER_FIELDS  

    
## -----------
class MIC_COADDListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=MIC_COADD  
    template_name = 'ddrug/mic_coadd/mic_coadd_list.html' 
    filterset_class=MIC_COADDfilter
    model_fields=model.HEADER_FIELDS

## -------------
class MIC_COADDCardView(MIC_COADDListView):
    template_name = 'ddrug/mic_coadd/mic_coadd_card.html'
  

## -----------
class MIC_PubListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=MIC_Pub  
    template_name = 'ddrug/mic_pub/mic_pub_list.html' 
    filterset_class=MIC_Pubfilter
    model_fields=model.HEADER_FIELDS

## -------------
class MIC_PubCardView(MIC_PubListView):
    template_name = 'ddrug/mic_pub/mic_pub_card.html'


## -------------
class BreakpointListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Breakpoint  
    template_name = 'ddrug/breakpoint/breakpoint_list.html' 
    filterset_class=Breakpointfilter
    model_fields=model.HEADER_FIELDS  

# --API Views--
## Drug
class API_Drug_List(API_ListView):
    queryset = Drug.objects.all()
    serializer_class = Drug_Serializer


## VITEK AST
class API_VITEK_ASTList(API_ListView):
    queryset = VITEK_AST.objects.all()
    serializer_class = VITEK_ASTSerializer

from rest_framework import generics
class API_Drug_Detail(generics.RetrieveAPIView):
    lookup_field = 'pk'
    queryset = Drug.objects.all()
    serializer_class = Drug_Serializer
 
# class VITEK_ASTCreate(generics.CreateAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer

# class VITEK_ASTUpdate(generics.RetrieveUpdateAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer

# class VITEK_ASTDelete(generics.DestroyAPIView):
#     queryset = VITEK_AST.objects.all()
#     serializer_class = VITEK_ASTSerializer
#
