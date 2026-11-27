from rest_framework import serializers
from rest_framework.exceptions import AuthenticationFailed
from rest_framework_simplejwt.serializers import TokenObtainPairSerializer

from .models import Drug, VITEK_AST


class AppUserTokenObtainPairSerializer(TokenObtainPairSerializer):
    """
    Block token issuance for users without application access.
    """

    def validate(self, attrs):
        data = super().validate(attrs)
        if not getattr(self.user, 'is_appuser', False):
            raise AuthenticationFailed(
                'User is not an application user',
                code='user_not_appuser',
            )
        return data


class Drug_Serializer(serializers.ModelSerializer):
    class Meta:
        model = Drug
        fields = '__all__'



class VITEK_ASTSerializer(serializers.ModelSerializer):
    class Meta:
        model = VITEK_AST
        fields = '__all__'