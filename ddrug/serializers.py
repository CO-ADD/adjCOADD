from rest_framework import serializers
from .models import Drug, VITEK_AST

class Drug_Serializer(serializers.ModelSerializer):
    class Meta:
        model = Drug
        fields = '__all__'



class VITEK_ASTSerializer(serializers.ModelSerializer):
    class Meta:
        model = VITEK_AST
        fields = '__all__'