import os
import json
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
from apputil.models import Dictionary, ApplicationUser
from apputil.utils.filters_base import FilteredListView
from apputil.utils.views_base import DeleteView_FKeyExist, ModelDeleteView
from apputil.utils.views_base import permission_not_granted
from .models import Gene
from .forms import (CreateGene_form, Genefilter,)

# --TAXONOMY Views--
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


@login_required
def createGene(req):
    kwargs={}
    kwargs['user']=req.user 
    form=CreateGene_form
    if req.method=='POST':
        form=CreateGene_form(req.POST)
        if form.is_valid():
            instance=form.save(commit=False)
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dgene/gene/gene_c.html', {'form':form})