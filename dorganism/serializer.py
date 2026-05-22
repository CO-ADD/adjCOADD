from django.contrib.auth.models import Group, User
from rest_framework import serializers

from dorganism.models import Organism, Organism_Batch

class Organism_Serializer(serializers.ModelSerializer):
    class Meta:
        model = Organism
        fields = ["organism_id", "organism_name", "pub_id",
                    'strain_ids','strain_code','strain_panel','strain_type',
                    'res_property','gen_property','seq_name','sero_clone','strain_identification',
                ]

class OrgBatch_Serializer(serializers.ModelSerializer):
    organism_id = Organism_Serializer()
    class Meta:
        model = Organism_Batch
        fields = ["orgbatch_id", "organism_id", 
                  "batch_id", "batch_notes", "batch_quality", "quality_source", "qc_status", "qc_record", "stock_date", "stock_level", 
                  ]
