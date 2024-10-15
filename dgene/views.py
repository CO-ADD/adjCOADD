import os
import json
from django_filters.views import FilterView

from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib import messages
from django.core.exceptions import ValidationError
from django.db import transaction, IntegrityError
from django.http import JsonResponse
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic.detail import DetailView
from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import FilteredListView
from apputil.utils.views_base import permission_not_granted, SimplecreateView, SimpleupdateView

from dgene.models import Genome_Sequence, ID_Pub, ID_Sequence, WGS_FastQC, WGS_CheckM, Gene, AMR_Genotype  
from dgene.forms import (GenomeSeq_Filter, GenomeSeq_Form,
                         WGS_FastQC_Filter, WGS_CheckM_Filter,
                         IDSeq_Filter, IDSeq_Form, IDPub_Form, IDPub_Filter,
                         Gene_Filter, Gene_Form, 
                         AMRGenotype_Filter,  
                         )



#=================================================================================================
# Genome Sequences
#=================================================================================================

class GenomeSeq_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= Genome_Sequence
    template_name = 'dgene/genomeseq/genomeseq_list.html' 
    filterset_class=GenomeSeq_Filter
    model_fields=model.HEADER_FIELDS

##
class GenomeSeq_CardView(GenomeSeq_ListView):
    template_name = 'dgene/genomeseq/genomeseq_card.html'

##
class GenomeSeq_CreateView(SimplecreateView):
    form_class=GenomeSeq_Form
    template_name='dgene/genomeseq/genomeseq_c.html'

##
class GenomeSeq_UpdateView(SimpleupdateView):
    form_class=GenomeSeq_Form
    template_name='dgene/genomeseq/genomeseq_u.html'
    model=ID_Sequence



#=================================================================================================
# ID Sequence 
#=================================================================================================

class IDSeq_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= ID_Sequence
    template_name = 'dgene/idseq/idseq_list.html' 
    filterset_class=IDSeq_Filter
    model_fields=model.HEADER_FIELDS

class IDPub_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= ID_Pub
    template_name = 'dgene/idseq/idpub_list.html' 
    filterset_class=IDPub_Filter
    model_fields=model.HEADER_FIELDS

##
class IDPub_CreateView(SimplecreateView):
    form_class=IDPub_Form
    template_name='dgene/idseq/idpub_c.html'

##
class IDPub_UpdateView(SimpleupdateView):
    form_class=IDPub_Form
    template_name='dgene/idseq/idpub_u.html'
    model=ID_Pub
    

#=================================================================================================
# WGS_FastQC - FastQ QC
#=================================================================================================
class WGS_FastQC_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= WGS_FastQC
    template_name = 'dgene/wgs_fastqc/fastqc_list.html' 
    filterset_class=WGS_FastQC_Filter
    model_fields=model.HEADER_FIELDS
    #ordering = []

#=================================================================================================
# WGS_CheckM - FastA CheckM
#=================================================================================================
class WGS_CheckM_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= WGS_CheckM
    template_name = 'dgene/wgs_checkm/checkm_list.html' 
    filterset_class=WGS_CheckM_Filter
    model_fields=model.HEADER_FIELDS
    #ordering = []


#=================================================================================================
# Genes
#=================================================================================================

class Gene_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= Gene
    template_name = 'dgene/gene/gene_list.html' 
    filterset_class=Gene_Filter
    model_fields=model.HEADER_FIELDS

##
class Gene_CardView(Gene_ListView):
    template_name = 'dgene/gene/gene_card.html'


class Gene_CreateView(SimplecreateView):
    form_class=Gene_Form
    template_name='dgene/gene/gene_c.html'

@login_required
def detailGene(req, pk):
    context={}
    object_=get_object_or_404(Gene, gene_id=pk)
    form=Gene_Form(instance=object_)    
    context["object"]=object_
    context["form"]=form
 
    return render(req, "dgene/gene/gene_detail.html", context)

class Gene_UpdateView(SimpleupdateView):
    form_class=Gene_Form
    template_name='dgene/gene/gene_u.html'
    model=Gene

#=================================================================================================
# AMR Genotype
#=================================================================================================

class AMRGenotype_ListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= AMR_Genotype
    template_name = 'dgene/gene/amrfinder_list.html' 
    filterset_class=AMRGenotype_Filter
    model_fields=model.HEADER_FIELDS
    ordering = ['orgbatch_id']
