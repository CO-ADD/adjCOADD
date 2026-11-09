from django.urls import path

from dbiologic.views import (Biologic_ListView, Biologic_CreateView, Biologic_DetailView, Biologic_UpdateView, Biologic_RemoveView,
                            Isoform_ListView, Isoform_CreateView, Isoform_DetailView, Isoform_UpdateView, Isoform_RemoveView,
                            BiologicModification_CreateView, BiologicModification_UpdateView,
                            BiologicModification_RemoveView, BiologicModification_ReorderView,
                            IsoBatch_CreateView, IsoBatch_UpdateView, IsoBatch_RemoveView,
                            Modification_ListView, Modification_CreateView, Modification_UpdateView,
                            Monomer_ListView, Monomer_CreateView, Monomer_UpdateView,
                            Sequence_StructureView)


urlpatterns = [
    # Biologic
    path('biologic_list', Biologic_ListView.as_view(), name="biologic_list"),
    path('biologic/<str:pk>', Biologic_DetailView, name="biologic_detail"),
    path('createBiologic/', Biologic_CreateView, name="biologic_create"),
    path('updateBiologic/<str:pk>', Biologic_UpdateView, name="biologic_update"),
    path('deleteBiologic/<str:pk>', Biologic_RemoveView.as_view(), name="biologic_delete"),

    # Biologic Isoform
    path('isoform_list', Isoform_ListView.as_view(), name="isoform_list"),
    path('isoform/<str:pk>', Isoform_DetailView, name="isoform_detail"),
    path('createIsoform/', Isoform_CreateView, name="isoform_create"),
    path('updateIsoform/<str:pk>', Isoform_UpdateView, name="isoform_update"),
    path('deleteIsoform/<str:pk>', Isoform_RemoveView.as_view(), name="isoform_delete"),

    # Isoform Modifications
    path('isoform/<str:pk>/modification/create', BiologicModification_CreateView, name="biologic_modification_create"),
    path('isoform/modification/update/<str:pk>', BiologicModification_UpdateView, name="biologic_modification_update"),
    path('isoform/modification/remove/<str:pk>', BiologicModification_RemoveView, name="biologic_modification_remove"),
    path('isoform/modification/reorder/<str:pk>/<str:direction>', BiologicModification_ReorderView, name="biologic_modification_reorder"),

    # Isoform Batch
    path('createIsoBatch/<str:isoform_id>', IsoBatch_CreateView, name="isobatch_create"),
    path('updateIsoBatch/<str:pk>', IsoBatch_UpdateView, name="isobatch_update"),
    path('deleteIsoBatch/<str:pk>', IsoBatch_RemoveView.as_view(), name="isobatch_delete"),

    # Modification Dictionary
    path('modification_list', Modification_ListView.as_view(), name="modification_list"),
    path('createModification/', Modification_CreateView, name="modification_create"),
    path('updateModification/<str:pk>', Modification_UpdateView, name="modification_update"),

    # Monomer
    path('monomer_list', Monomer_ListView.as_view(), name="monomer_list"),
    path('createMonomer/', Monomer_CreateView, name="monomer_create"),
    path('updateMonomer/<str:pk>', Monomer_UpdateView, name="monomer_update"),

    # Sequence structure preview
    path('sequence_structure', Sequence_StructureView, name="sequence_structure"),
]
