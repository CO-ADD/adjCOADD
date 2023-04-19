from django.conf import settings
from django.conf.urls.static import static
from django.urls import path, include, re_path

from .views import  (GeneCardView,GeneListView,detailGene,
                     createGene,updateGene)
urlpatterns = [
    # Taxonomy
    path('gene_card', GeneCardView.as_view(), name='gene_card'),
    path('gene_list', GeneListView.as_view(), name="gene_list"),
    path('gene/<str:pk>', detailGene, name="gene_detail"),
    path('createGene/', createGene, name="gene_create"),
    path('updateGene/<str:pk>', updateGene, name="gene_update"),
    # path('deleteTax/<slug:slug>', deleteTaxonomy, name="taxonomy_delete"),
]