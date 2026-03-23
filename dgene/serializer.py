
from django.contrib.auth.models import Group, User
from rest_framework import serializers

from dgene.models import Genome_Sequence, ID_Sequence
from dorganism.serializer import OrgBatch_Serializer

class GenomeSeq_Serializer(serializers.ModelSerializer):
    
    #tag_detail = serializers.HyperlinkedIdentityField(view_name='genomeseq_api_list')
    orgbatch_id = OrgBatch_Serializer()
    class Meta:
        model = Genome_Sequence
        fields = ["seq_id", "run_id", "orgbatch_id", "seq_code",
                  "seq_type","seq_method","seq_files",
                  "runsample_file","runsample_dir","runsample_name",
                #  "tag_detail",
                  ]
        
class IDSeq_Serializer(serializers.ModelSerializer):
    
    class Meta:
        model = ID_Sequence
        fields = ["seq_id", "seq_file","seq_filetype", "kraken_organisms",
                "mlst_scheme","mlst_seqtype","mlst_alleles",
                "gtdbtk_class","gtdbtk_fastani",
                "source","id_notes",
        ]
