from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (GeneCardView,GeneListView,detailGene,GeneCreateView,GeneUpdateView,
                     ID_PubListView,ID_PubCreateView,ID_PubUpdateView,
                     SequenceCardView,SequenceListView,SequenceCreateView,SequenceUpdateView, 
                     WGS_FastQCListView, WGS_CheckMListView)
                     
urlpatterns = [
    # gene
    path('gene_list', GeneListView.as_view(), name="gene_list"),
    path('gene/<str:pk>', detailGene, name="gene_detail"),
    path('createGene/', GeneCreateView.as_view(), name="gene_create"),
    path('updateGene/<str:pk>', GeneUpdateView.as_view(), name="gene_update"),

    # id_pub
    path('id_pub_list', ID_PubListView.as_view(), name="id_pub_list"),
    # path('id_pub/<str:pk>', detailGene, name="id_pub_detail"),
    path('createid_pub/', ID_PubCreateView.as_view(), name="id_pub_create"),
    path('updateid_pub/<str:pk>', ID_PubUpdateView.as_view(), name="id_pub_update"),

    # sequence
    path('sequence_list', SequenceListView.as_view(), name="sequence_list"),
    path('createSequence/', SequenceCreateView.as_view(), name="sequence_create"),
    path('updateSequence/<str:pk>', SequenceUpdateView.as_view(), name="sequence_update"),

    #wgs
    path('wgs-fastqc_list', WGS_FastQCListView.as_view(), name="wgs-fastqc_list"),
    path('wgs-checkm_list', WGS_CheckMListView.as_view(), name="wgs-checkm_list"),


]