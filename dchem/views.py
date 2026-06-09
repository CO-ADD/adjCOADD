from django.shortcuts import render

from rest_framework import viewsets, filters
from django_rdkit.models import MOL_FROM_SMILES, TANIMOTO_SMILES

from dchem.models import Chem_Structure
from dchem.serializer import ChemStructure_Serializer

# -----------------------------------------------------------------
class Compound_ListAPI(viewsets.ModelViewSet):
# -----------------------------------------------------------------
    queryset = Chem_Structure.objects.all()
    serializer_class = ChemStructure_Serializer

    def get_queryset(self):
        queryset = Chem_Structure.objects.all()
        
        # Substructure Filter
        substructure = self.request.query_params.get('substruct')
        if substructure:
            q_mol = MOL_FROM_SMILES(substructure)
            queryset = queryset.filter(molecule__contains=q_mol)

        # Similarity Filter
        similarity_smile = self.request.query_params.get('similar_to')
        threshold = float(self.request.query_params.get('threshold', 0.85))
        if similarity_smile:
            queryset = queryset.annotate(
                similarity=TANIMOTO_SMILES('molecule', similarity_smile)
            ).filter(similarity__gte=threshold).order_by('-similarity')

        return queryset
    

