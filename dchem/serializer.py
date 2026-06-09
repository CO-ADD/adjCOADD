from rest_framework import serializers
from dchem.models import Chem_Structure

from rdkit import Chem

# -----------------------------------------------------------------
class ChemStructure_Serializer(serializers.ModelSerializer):
# -----------------------------------------------------------------
    smiles = serializers.SerializerMethodField()
    molblock = serializers.SerializerMethodField()

    class Meta:
        model = Chem_Structure
        fields = ['structure_id', 'smiles', 'molblock']

    def get_smiles(self, obj):
        return Chem.MolToSmiles(obj.molecule)

    def get_molblock(self, obj):
        return Chem.MolToMolBlock(obj.molecule)
