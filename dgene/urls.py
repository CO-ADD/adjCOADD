from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from dgene.views import  (GenomeSeq_ListView,GenomeSeq_CreateView,GenomeSeq_UpdateView, 
                     IDSeq_ListView,
                     WGS_FastQC_ListView, WGS_CheckM_ListView,
                     Gene_ListView,Gene_CardView,detailGene,Gene_CreateView,Gene_UpdateView,
                     AMRGenotype_ListView,
                     IDPub_ListView,IDPub_CreateView,IDPub_UpdateView,
                     )
                     
urlpatterns = [
    # sequence
    path('sequence_list', GenomeSeq_ListView.as_view(), name="genomeseq_list"),
    path('sequence/<str:pk>', detailGene, name="gene_detail"),
    path('createSequence/', GenomeSeq_CreateView.as_view(), name="genomeseq_create"),
    #path('updateSequence/<str:pk>', GenomeSeq_UpdateView.as_view(), name="genomeseq_update"),

    #wgs
    path('wgs_fastqc_list', WGS_FastQC_ListView.as_view(), name="wgs_fastqc_list"),
    path('wgs_checkm_list', WGS_CheckM_ListView.as_view(), name="wgs_checkm_list"),

    # ID by Seq
    path('idseq_list', IDSeq_ListView.as_view(), name="idseq_list"),
    #path('createSequence/', SequenceCreateView.as_view(), name="sequence_create"),
    #path('updateSequence/<str:pk>', SequenceUpdateView.as_view(), name="sequence_update"),

    # gene
    path('gene_list', Gene_ListView.as_view(), name="gene_list"),
    path('gene/<str:pk>', detailGene, name="gene_detail"),
    path('createGene/', Gene_CreateView.as_view(), name="gene_create"),
    #path('updateGene/<str:pk>', GeneUpdateView.as_view(), name="gene_update"),

    # gene
    path('amrgene_list', AMRGenotype_ListView.as_view(), name="amrgene_list"),

    # id_pub
    path('idpub_list', IDPub_ListView.as_view(), name="id_pub_list"),
    # path('id_pub/<str:pk>', detailGene, name="id_pub_detail"),
    path('createid_pub/', IDPub_CreateView.as_view(), name="id_pub_create"),
    path('updateid_pub/<str:pk>', IDPub_UpdateView.as_view(), name="id_pub_update"),

]