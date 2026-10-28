rom django.urls import path

from dchem.views import (Compound_ListAPI, Chem_ListView, Chem_DetailView, Chem_CreateView, Chem_UpdateView)


urlpatterns = [
    # Chem Structure API
    path('api/chemstructure', Compound_ListAPI.as_view({'get': 'list'}), name="chemstructure_api_list"),

    # Chem Structure (filter-first, large table)
    path('chemstructure_list', Chem_ListView.as_view(key='structure'), name="chemstructure_list"),
    path('structure/<str:pk>', Chem_DetailView, {'key': 'structure'}, name="chemstructure_detail"),

    # Chem Salt
    path('chemsalt_list', Chem_ListView.as_view(key='salt'), name="chemsalt_list"),
    path('createSalt/', Chem_CreateView, {'key': 'salt'}, name="chemsalt_create"),
    path('updateSalt/<str:pk>', Chem_UpdateView, {'key': 'salt'}, name="chemsalt_update"),

    # Chem Group
    path('chemgroup_list', Chem_ListView.as_view(key='group'), name="chemgroup_list"),
    path('createChemGroup/', Chem_CreateView, {'key': 'group'}, name="chemgroup_create"),
    path('updateChemGroup/<str:pk>', Chem_UpdateView, {'key': 'group'}, name="chemgroup_update"),

    # Chem Reaction
    path('chemreaction_list', Chem_ListView.as_view(key='reaction'), name="chemreaction_list"),
    path('createChemReaction/', Chem_CreateView, {'key': 'reaction'}, name="chemreaction_create"),
    path('updateChemReaction/<str:pk>', Chem_UpdateView, {'key': 'reaction'}, name="chemreaction_update"),

    # Amino Acid
    path('aminoacid_list', Chem_ListView.as_view(key='aminoacid'), name="aminoacid_list"),
    path('aminoacid/<str:pk>', Chem_DetailView, {'key': 'aminoacid'}, name="aminoacid_detail"),
    path('createAminoAcid/', Chem_CreateView, {'key': 'aminoacid'}, name="aminoacid_create"),
    path('updateAminoAcid/<str:pk>', Chem_UpdateView, {'key': 'aminoacid'}, name="aminoacid_update"),
]
