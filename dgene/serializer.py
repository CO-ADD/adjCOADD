
from django.contrib.auth.models import Group, User
from rest_framework import serializers

from dgene.models import Genome_Sequence

class GenomeSeq_Serializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Genome_Sequence
        fields = ["seq_id", "run_id", "orgbatch_id", "seq_type","seq_method","seq_files",
                  "runsample_file","runsample_dir","runsample_name"]
        
