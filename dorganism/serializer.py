from django.contrib.auth.models import Group, User
from rest_framework import serializers

from dorganism.models import Organism, Organism_Batch

class Organism_Serializer(serializers.ModelSerializer):
    
    class Meta:
        model = Organism
        fields = ["organism_id", "organism_name", "pub_id",
                  ]

class OrgBatch_Serializer(serializers.ModelSerializer):
    organism_id = Organism_Serializer()
    class Meta:
        model = Organism_Batch
        fields = ["orgbatch_id", "organism_id"
                  ]
