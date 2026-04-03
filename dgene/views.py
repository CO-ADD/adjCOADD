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

from rest_framework import viewsets, status
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from django_filters.rest_framework import DjangoFilterBackend

from apputil.models import Dictionary, ApplicationUser
from applib.django.base.views import permission_not_granted, Base_CreateView, Base_UpdateView, Filtered_ListView

from dgene.models import Genome_Sequence, ID_Pub, ID_Sequence, WGS_FastQC, WGS_CheckM, Gene, AMR_Genotype  
from dgene.forms import (GenomeSeq_Filter, GenomeSeq_Form,
                         WGS_FastQC_Filter, WGS_CheckM_Filter,
                         IDSeq_Filter, IDSeq_Form, IDPub_Form, IDPub_Filter,
                         Gene_Filter, Gene_Form, 
                         AMRGenotype_Filter,  
                         )
from dgene.serializer import GenomeSeq_Serializer, IDSeq_Serializer, WGS_CheckM_Serializer



#=================================================================================================
# Genome Sequences
#=================================================================================================

class GenomeSeq_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= Genome_Sequence
    template_name = 'dgene/genomeseq/genomeseq_list.html' 
    filterset_class=GenomeSeq_Filter
    model_fields=model.LIST_VIEW_FIELDS
    ordering = ['seq_id']

##
class GenomeSeq_CardView(GenomeSeq_ListView):
    template_name = 'dgene/genomeseq/genomeseq_card.html'

##
class GenomeSeq_CreateView(Base_CreateView):
    form_class=GenomeSeq_Form
    template_name='dgene/genomeseq/genomeseq_c.html'

##
class GenomeSeq_UpdateView(Base_UpdateView):
    form_class=GenomeSeq_Form
    template_name='dgene/genomeseq/genomeseq_u.html'
    model=ID_Sequence

##
class GenomeSeq_ListAPI(viewsets.ModelViewSet):
    queryset = Genome_Sequence.objects.all()
    serializer_class = GenomeSeq_Serializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['seq_id','run_id','orgbatch_id','seq_code',]


class GenomeSeq_DetailAPI(viewsets.ModelViewSet):
    queryset = Genome_Sequence.objects.all()
    serializer_class = GenomeSeq_Serializer

    permission_classes = [IsAuthenticated]
    
    # def get_queryset(self):
    #     obj= Genome_Sequence.objects.get(seq_id=self.kwargs['pk'])
    #     user = self.request.user
    #     print(f" [get_queryset] {obj}")
    #     return super().get_queryset()

    def update(self, request, *args, **kwargs):
        obj= Genome_Sequence.objects.get(seq_id=self.kwargs['pk'])
        user = self.request.user
        serializer = self.serializer_class(obj,data=request.data, partial=True)
        if serializer.is_valid():
            serializer.save()
            print(f" [update] Genome_Sequence {obj} <- {request.data} {serializer.is_valid()}")
            updated_data = serializer.data
            updated_data['custom_message'] = f"{str(obj)} updated successfully!"
            return Response(updated_data, status=status.HTTP_200_OK)
        Response(serializer.errors, status=400)

#=================================================================================================
# ID Sequence 
#=================================================================================================

class IDSeq_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= ID_Sequence
    template_name = 'dgene/idseq/idseq_list.html' 
    filterset_class=IDSeq_Filter
    model_fields=model.LIST_VIEW_FIELDS

##
class IDSeq_ListAPI(viewsets.ModelViewSet):
    queryset = ID_Sequence.objects.all()
    serializer_class = IDSeq_Serializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['seq_id','seq_file']


class IDSeq_UpdateAPI(viewsets.ModelViewSet):
    queryset = ID_Sequence.objects.all()
    serializer_class = IDSeq_Serializer

    permission_classes = [IsAuthenticated]
    
    # def get_queryset(self):
    #     obj= Genome_Sequence.objects.get(seq_id=self.kwargs['pk'])
    #     user = self.request.user
    #     print(f" [get_queryset] {obj}")
    #     return super().get_queryset()

    def update(self, request, *args, **kwargs):
        
        seq_id = request.data.get('seq_id',None)
        seq_filetype = request.data.get('seq_filetype',None)
        seq_file = request.data.get('seq_file',None)
        
        qryDict = {
            'seq_id' : seq_id,
            'seq_filetype':seq_filetype,
            'seq_file':seq_file
        }
        
        obj=ID_Sequence.get(seq_id,seq_filetype,seq_file)
        if obj:
            serializer = self.serializer_class(obj,data=request.data, partial=True)
            if serializer.is_valid():
                serializer.save()
                #print(f" [update] IDSeq_UpdateAPI {obj} <- {request.data} {serializer.is_valid()}")
                updated_data = serializer.data
                updated_data['custom_message'] = f"{str(obj)} updated successfully!"
                return Response(updated_data, status=status.HTTP_200_OK)
        else: 
            serializer = self.get_serializer(data=request.data)
            if serializer.is_valid():
                self.perform_create(serializer)
                updated_data = serializer.data
                updated_data['custom_message'] = f"{str(obj)} created successfully!"
                return Response(updated_data, status=status.HTTP_201_CREATED)
            else:
                print(serializer.errors)
        Response(serializer.errors, status=400)

#=================================================================================================
# ID Public
#=================================================================================================
class IDPub_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= ID_Pub
    template_name = 'dgene/idpub/idpub_list.html' 
    filterset_class=IDPub_Filter
    model_fields=model.LIST_VIEW_FIELDS

##
class IDPub_CreateView(Base_CreateView):
    form_class=IDPub_Form
    template_name='dgene/idpub/idpub_c.html'

##
class IDPub_UpdateView(Base_UpdateView):
    form_class=IDPub_Form
    template_name='dgene/idpub/idpub_u.html'
    model=ID_Pub
    

#=================================================================================================
# WGS_FastQC - FastQ QC
#=================================================================================================
class WGS_FastQC_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= WGS_FastQC
    template_name = 'dgene/wgs_fastqc/fastqc_list.html' 
    filterset_class=WGS_FastQC_Filter
    model_fields=model.LIST_VIEW_FIELDS
    #ordering = []

#=================================================================================================
# WGS_CheckM - FastA CheckM
#=================================================================================================
class WGS_CheckM_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= WGS_CheckM
    template_name = 'dgene/wgs_checkm/checkm_list.html' 
    filterset_class=WGS_CheckM_Filter
    model_fields=model.LIST_VIEW_FIELDS
    #ordering = []

class WGS_CheckM_ListAPI(viewsets.ModelViewSet):
    queryset = WGS_CheckM.objects.all()
    serializer_class = WGS_CheckM_Serializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['seq_id','assembler','fasta_type']

class WGS_CheckM_UpdateAPI(viewsets.ModelViewSet):
    queryset = WGS_CheckM.objects.all()
    serializer_class = WGS_CheckM_Serializer

    permission_classes = [IsAuthenticated]
    
    print(f" [WGS_CheckM_UpdateAPI]")
    # def get_queryset(self):
    #     obj= Genome_Sequence.objects.get(seq_id=self.kwargs['pk'])
    #     user = self.request.user
    #     print(f" [get_queryset] {obj}")
    #     return super().get_queryset()

    def update(self, request, *args, **kwargs):
        
        seq_id = request.data.get('seq_id',None)
        assembler = request.data.get('assembler',None)
        fasta_type = request.data.get('fasta_type',None)
        
        qryDict = {
            'seq_id' : seq_id,
            'assembler':assembler,
            'fasta_type':fasta_type
        }
        
        print(f" [WGS_CheckM_UpdateAPI] {request.data}")
        obj=WGS_CheckM.get(seq_id,assembler,fasta_type)
        print(f" [WGS_CheckM_UpdateAPI] {obj}")
        if obj:
            serializer = self.serializer_class(obj,data=request.data, partial=True)
            if serializer.is_valid():
                serializer.save()
                print(f" [update] WGS_CheckM_UpdateAPI {obj} <- {request.data} {serializer.is_valid()}")
                updated_data = serializer.data
                updated_data['custom_message'] = f"{str(obj)} updated successfully!"
                return Response(updated_data, status=status.HTTP_200_OK)
        else: 
            serializer = self.get_serializer(data=request.data)
            if serializer.is_valid():
                print(f" [create] WGS_CheckM_UpdateAPI <- {request.data} {serializer.is_valid()}")
                self.perform_create(serializer)
                updated_data = serializer.data
                updated_data['custom_message'] = f"{str(obj)} created successfully!"
                return Response(updated_data, status=status.HTTP_201_CREATED)
            else:
                print(serializer.errors)
        Response(serializer.errors, status=400)

#=================================================================================================
# Genes
#=================================================================================================

class Gene_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= Gene
    template_name = 'dgene/gene/gene_list.html' 
    filterset_class=Gene_Filter
    model_fields=model.LIST_VIEW_FIELDS

##
class Gene_CardView(Gene_ListView):
    template_name = 'dgene/gene/gene_card.html'

class Gene_CreateView(Base_CreateView):
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

class Gene_UpdateView(Base_UpdateView):
    form_class=Gene_Form
    template_name='dgene/gene/gene_u.html'
    model=Gene

#=================================================================================================
# AMR Genotype
#=================================================================================================

class AMRGenotype_ListView(LoginRequiredMixin, Filtered_ListView):
    login_url = '/'
    model= AMR_Genotype
    template_name = 'dgene/amrgenotype/amrgenotype_list.html' 
    filterset_class=AMRGenotype_Filter
    model_fields=model.LIST_VIEW_FIELDS
    ordering = ['orgbatch_id']
