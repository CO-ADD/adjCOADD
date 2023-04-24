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
from .models import Gene, ID_Pub, ID_Sequence, WGS_FastQC, WGS_CheckM
from .forms import (Gene_form, Genefilter, Sequence_form, Sequencefilter, FastQCfilter, CheckMfilter, ID_Pub_form, ID_Pubfilter)

# --Gene Views--
##
class GeneListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= Gene
    template_name = 'dgene/gene/gene_list.html' 
    filterset_class=Genefilter
    model_fields=model.HEADER_FIELDS

##
class GeneCardView(GeneListView):
    template_name = 'dgene/gene/gene_card.html'


class GeneCreateView(SimplecreateView):
    form_class=Gene_form
    template_name='dgene/gene/gene_c.html'

@login_required
def detailGene(req, pk):
    context={}
    object_=get_object_or_404(Gene, gene_id=pk)
    form=Gene_form(instance=object_)    
    context["object"]=object_
    context["form"]=form
 
    return render(req, "dgene/gene/gene_detail.html", context)


class GeneUpdateView(SimpleupdateView):
    form_class=Gene_form
    template_name='dgene/gene/gene_u.html'
    model=Gene

# -- ID-Pub Views--
##
class ID_PubListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= ID_Pub
    template_name = 'dgene/id_pub/id_pub_list.html' 
    filterset_class=ID_Pubfilter
    model_fields=model.HEADER_FIELDS

##
class ID_PubCreateView(SimplecreateView):
    form_class=ID_Pub_form
    template_name='dgene/id_pub/id_pub_c.html'

##
class ID_PubUpdateView(SimpleupdateView):
    form_class=ID_Pub_form
    template_name='dgene/id_pub/id_pub_u.html'
    model=ID_Sequence
    
# --Sequence Views--
##
class SequenceListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= ID_Sequence
    template_name = 'dgene/id_sequence/sequence_list.html' 
    filterset_class=Sequencefilter
    model_fields=model.HEADER_FIELDS

##
class SequenceCardView(SequenceListView):
    template_name = 'dgene/id_sequence/sequence_card.html'

##
class SequenceCreateView(SimplecreateView):
    form_class=Sequence_form
    template_name='dgene/id_sequence/sequence_c.html'

##
class SequenceUpdateView(SimpleupdateView):
    form_class=Sequence_form
    template_name='dgene/id_sequence/sequence_u.html'
    model=ID_Sequence

# --WGS FastQC--
##
class WGS_FastQCListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= WGS_FastQC
    template_name = 'dgene/wgs_fastqc/fastqc_list.html' 
    filterset_class=FastQCfilter
    model_fields=model.HEADER_FIELDS

# --WGS CheckM--
##
class WGS_CheckMListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model= WGS_CheckM
    template_name = 'dgene/wgs_checkm/checkm_list.html' 
    filterset_class=CheckMfilter
    model_fields=model.HEADER_FIELDS

