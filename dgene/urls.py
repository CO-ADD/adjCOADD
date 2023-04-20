from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (GeneCardView,GeneListView,detailGene,GeneCreateView,GeneUpdateView,
                     SequenceCardView,SequenceListView,SequenceCreateView,SequenceUpdateView, 
                     WGS_FastQCListView, WGS_FastQCCardView, WGS_CheckMListView, WGS_CheckMCardView)
                     
urlpatterns = [
    # gene
    path('gene_card', GeneCardView.as_view(), name='gene_card'),
    path('gene_list', GeneListView.as_view(), name="gene_list"),
    path('gene/<str:pk>', detailGene, name="gene_detail"),
    path('createGene/', GeneCreateView.as_view(), name="gene_create"),
    path('updateGene/<str:pk>', GeneUpdateView.as_view(), name="gene_update"),

    # sequence
    path('sequence_card', SequenceCardView.as_view(), name='sequence_card'),
    path('sequence_list', SequenceListView.as_view(), name="sequence_list"),
    path('createSequence/', SequenceCreateView.as_view(), name="sequence_create"),
    path('updateSequence/<str:pk>', SequenceUpdateView.as_view(), name="sequence_update"),

    #wgs
    path('wgs-fastqc_card', WGS_FastQCCardView.as_view(), name='wgs-fastqc_card'),
    path('wgs-fastqc_list', WGS_FastQCListView.as_view(), name="wgs-fastqc_list"),
    path('wgs-checkm_card', WGS_CheckMCardView.as_view(), name='wgs-checkm_card'),
    path('wgs-checkm_list', WGS_CheckMListView.as_view(), name="wgs-checkm_list"),


]