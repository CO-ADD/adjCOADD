from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dgene.views import  (GenomeSeq_ListView,GenomeSeq_CreateView,GenomeSeq_UpdateView, GenomeSeq_ListAPI, GenomeSeq_DetailAPI,
                     IDSeq_ListView, IDSeq_ListAPI, IDSeq_UpdateAPI,
                     WGS_FastQC_ListView, 
                     WGS_CheckM_ListView, WGS_CheckM_ListAPI, WGS_CheckM_UpdateAPI,
                     Gene_ListView,Gene_CardView,detailGene,Gene_CreateView,Gene_UpdateView, Gene_ListAPI, Gene_UpdateAPI,
                     AMRGenotype_ListView, AMRGenotype_ListAPI, AMRGenotype_UpdateAPI,
                     IDPub_ListView,IDPub_CreateView,IDPub_UpdateView,
                     )
                     
urlpatterns = [

    # Genome_Sequence
    path('sequence_list', GenomeSeq_ListView.as_view(), name="genomeseq_list"),
    path('sequence/<str:pk>', detailGene, name="gene_detail"),
    path('createSequence/', GenomeSeq_CreateView.as_view(), name="genomeseq_create"),
    #path('updateSequence/<str:pk>', GenomeSeq_UpdateView.as_view(), name="genomeseq_update"),
    path('api/sequence_list', GenomeSeq_ListAPI.as_view({'get': 'list'}), name="genomeseq_api_list"),
    path('api/sequence/<str:pk>', GenomeSeq_DetailAPI.as_view({'get': 'retrieve',"patch": "partial_update","post": "update"}), name="genomeseq_api_detail"),
    

    # WGS_FastQC
    path('wgs_fastqc_list', WGS_FastQC_ListView.as_view(), name="wgs_fastqc_list"),
 
    # WGS_CheckM
    path('wgs_checkm_list', WGS_CheckM_ListView.as_view(), name="wgs_checkm_list"),
    path('api/checkm_list', WGS_CheckM_ListAPI.as_view({'get': 'list'}), name="checkm_api_list"),
    path('api/checkm_update', WGS_CheckM_UpdateAPI.as_view({"patch": "partial_update","post": "update"}), name="checkm_api_update"),
 
    # ID_Sequence
    path('idseq_list', IDSeq_ListView.as_view(), name="idseq_list"),
    path('api/idseq_list', IDSeq_ListAPI.as_view({'get': 'list'}), name="idseq_api_list"),
    path('api/idseq_update', IDSeq_UpdateAPI.as_view({"patch": "partial_update","post": "update"}), name="idseq_api_detail"),
    #path('createSequence/', SequenceCreateView.as_view(), name="sequence_create"),
    #path('updateSequence/<str:pk>', SequenceUpdateView.as_view(), name="sequence_update"),

    # ID_Pub
    path('idpub_list', IDPub_ListView.as_view(), name="id_pub_list"),
    # path('id_pub/<str:pk>', detailGene, name="id_pub_detail"),
    path('createid_pub/', IDPub_CreateView.as_view(), name="id_pub_create"),
    path('updateid_pub/<str:pk>', IDPub_UpdateView.as_view(), name="id_pub_update"),

    # Gene
    path('gene_list', Gene_ListView.as_view(), name="gene_list"),
    path('gene/<str:pk>', detailGene, name="gene_detail"),
    path('createGene/', Gene_CreateView.as_view(), name="gene_create"),
    path('api/gene_list', Gene_ListAPI.as_view({'get': 'list'}), name="gene_api_list"),
    path('api/gene_update', Gene_UpdateAPI.as_view({"patch": "partial_update","post": "update"}), name="gene_api_detail"),
    #path('updateGene/<str:pk>', GeneUpdateView.as_view(), name="gene_update"),

    # AMR_Genotype
    path('amrgene_list', AMRGenotype_ListView.as_view(), name="amrgene_list"),
    path('api/amrgene_list', AMRGenotype_ListAPI.as_view({'get': 'list'}), name="idseq_api_list"),
    path('api/amrgene_update', AMRGenotype_UpdateAPI.as_view({"patch": "partial_update","post": "update"}), name="idseq_api_detail"),


]