rom django.urls import path

from dpeptide.views import (Peptide_ListView, Peptide_CreateView, Peptide_DetailView, Peptide_UpdateView, Peptide_RemoveView,
                            Isoform_ListView, Isoform_CreateView, Isoform_DetailView, Isoform_UpdateView, Isoform_RemoveView,
                            IsoBatch_CreateView, IsoBatch_UpdateView, IsoBatch_RemoveView,
                            Modification_ListView, Modification_CreateView, Modification_UpdateView)


urlpatterns = [
    # Peptide
    path('peptide_list', Peptide_ListView.as_view(), name="peptide_list"),
    path('peptide/<str:pk>', Peptide_DetailView, name="peptide_detail"),
    path('createPeptide/', Peptide_CreateView, name="peptide_create"),
    path('updatePeptide/<str:pk>', Peptide_UpdateView, name="peptide_update"),
    path('deletePeptide/<str:pk>', Peptide_RemoveView.as_view(), name="peptide_delete"),

    # Peptide Isoform
    path('isoform_list', Isoform_ListView.as_view(), name="isoform_list"),
    path('isoform/<str:pk>', Isoform_DetailView, name="isoform_detail"),
    path('createIsoform/', Isoform_CreateView, name="isoform_create"),
    path('updateIsoform/<str:pk>', Isoform_UpdateView, name="isoform_update"),
    path('deleteIsoform/<str:pk>', Isoform_RemoveView.as_view(), name="isoform_delete"),

    # Isoform Batch
    path('createIsoBatch/<str:isoform_id>', IsoBatch_CreateView, name="isobatch_create"),
    path('updateIsoBatch/<str:pk>', IsoBatch_UpdateView, name="isobatch_update"),
    path('deleteIsoBatch/<str:pk>', IsoBatch_RemoveView.as_view(), name="isobatch_delete"),

    # Modification Dictionary
    path('modification_list', Modification_ListView.as_view(), name="modification_list"),
    path('createModification/', Modification_CreateView, name="modification_create"),
    path('updateModification/<str:pk>', Modification_UpdateView, name="modification_update"),
]
