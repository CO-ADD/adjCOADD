from rest_framework import serializers
from .models import VITEK_AST

class VITEK_ASTSerializer(serializers.ModelSerializer):
    class Meta:
        model = VITEK_AST
        fields = '__all__'