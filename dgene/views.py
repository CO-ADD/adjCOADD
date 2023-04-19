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
from apputil.utils.views_base import DeleteView_FKeyExist, ModelDeleteView
from apputil.utils.views_base import permission_not_granted
from .models import Gene
from .forms import (Gene_form, Genefilter)

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
    form=Gene_form
    if req.method=='POST':
        form=Gene_form(req.POST)
        if form.is_valid():
            instance=form.save(commit=False)
            instance.save(**kwargs)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            messages.error(req, form.errors)
            return redirect(req.META['HTTP_REFERER'])      
    return render(req, 'dgene/gene/gene_c.html', {'form':form})

@login_required
def detailGene(req, pk):
    context={}
    object_=get_object_or_404(Gene, gene_id=pk)
    form=Gene_form(instance=object_)
    # Array display handler
    # if object_.organisms:
    #     context["gene_organisms"]=",".join(object_.organisms)
    # else:
    #     context["gene_organisms"]=""
    
    context["object"]=object_
    context["form"]=form
 
    return render(req, "dgene/gene/gene_detail.html", context)


@login_required
def updateGene(req, pk):
    object_=get_object_or_404(Gene, pk=pk)
    kwargs={}
    kwargs['user']=req.user 
    form=Gene_form(instance=object_)
    if req.method=='POST':
        form=Gene_form(req.POST, instance=object_)
        if form.is_valid():
            instance=form.save(commit=False)
            try:        
                instance.save(**kwargs)
            except Exception as err:
                messages.error(req, err)
            return redirect(req.META['HTTP_REFERER']) 
        else:
            print(form.errors)
    return render(req, 'ddrug/drug/drug_u.html', {'form':form, 'object':object_})