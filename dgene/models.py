from django.db import models

from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary
from dorganism.models import Organism, Organism_Batch
from dscreen.models import Screen_Run

#-------------------------------------------------------------------------------------------------
# Gene and Identificationrelated Application Model
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Gene(AuditModel):
    """
    List of Genes
    """
#=================================================================================================
    HEADER_FIELDS = {
        "gene_name":"Gene Name",
        "gene_othernames":"Other Name",
        "gene_type":"Gene Type",
        "gene_class":"Gene Class",
        "gene_note":"Gene Note",
    }

    Choice_Dictionary = {
        'gene_type':'Gene_Type',
    }

    ID_SEQUENCE = 'Gene'
    ID_PREFIC = 'GEN'
    ID_PAD = 5

    gene_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Drug ID")
    gene_name = models.CharField(max_length=25, unique=True,  verbose_name = "Gene Name")
    urlname = models.SlugField(max_length=30, verbose_name = "URLGene")
    gene_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Gene Type", on_delete=models.DO_NOTHING,
         db_column="gene_type", related_name="%(class)s_genetype")
    gene_note = models.CharField(max_length=120, verbose_name = "Gene Note")
    gene_modification = models.CharField(max_length=120, verbose_name = "Modification")
    protein_name = models.CharField(max_length=150, blank=True,  verbose_name = "Protein")
    protein_class = models.CharField(max_length=50, blank=True,  verbose_name = "Class")
    protein_subclass = models.CharField(max_length=50,  blank=True, verbose_name = "SubClass")
    resistance_class = models.CharField(max_length=50, blank=True,  verbose_name = "Resistance Class")
    variants = models.CharField(max_length=120, blank=True,  verbose_name = "Variants")
    genbank = models.CharField(max_length=120, blank=True,  verbose_name = "GenBank")
    uniprot = models.CharField(max_length=120, blank=True,  verbose_name = "UniProt")
    organisms =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=10, verbose_name = "Organisms", null=True, blank=True)
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'gene'
        ordering=['gene_type','resistance_class','gene_name']
        indexes = [
             models.Index(name="gene_gtype_idx",fields=['gene_type']),
             models.Index(name="gene_rclass_idx",fields=['resistance_class']),
             models.Index(name="gene_gmod_idx",fields=['gene_modification']),
             models.Index(name="gene_pclass_idx",fields=['protein_class']),
             models.Index(name="gene_psclass_idx",fields=['protein_subclass']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.gene_name}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.gene_name} {self.resistance_class} {self.gene_type}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,GeneName,verbose=0):
    # Returns an instance if found by [GeneName]
        try:
            retInstance = cls.objects.get(gene_name=GeneName)
        except:
            if verbose:
                print(f"[Gene Not Found] {GeneName} ")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,GeneName,verbose=0):
    # Returns an instance if found by [GeneName]
        return cls.objects.filter(gene_name=GeneName).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.gene_id:
            self.gene_id = self.next_id()
            if self.gene_id: 
                super(Gene, self).save(*args, **kwargs)
        else:
            super(Gene, self).save(*args, **kwargs) 


#=================================================================================================
class ID_Pub(AuditModel):
    """
     Identification from Public or Collaborative sources    
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "organism_id":"Organism ID",
        "id_type":"ID Method",
        "id_organism":"Organism",
        #"id_probability": "Probability",
        #"id_confidence": "Confidence",
        "source": "Source",
        "id_date":"Date",
    }
    Choice_Dictionary = {
        'id_type':'ID_Type',
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id") 
    id_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "ID Method", on_delete=models.DO_NOTHING,
         db_column="id_type", related_name="%(class)s_idtype")
    id_organism = models.CharField(max_length=120,  blank=True, verbose_name = "ID Organism")
    #id_probability = models.CharField(max_length=120, blank=True,  verbose_name = "ID Probability")
    #id_confidence = models.CharField(max_length=120, blank=True,  verbose_name = "ID Confidence")
    id_date = models.DateField(blank=True, verbose_name = "ID Date")
    source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'id_pub'
        ordering=['organism_id','id_type']
        indexes = [
             models.Index(name="idp_org_idx",fields=['id_organism']),
             models.Index(name="idp_idtype_idx",fields=['id_type']),
             models.Index(name="idp_source_idx",fields=['source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.organism_id} {self.id_organism} "
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.organism_id} {self.id_organism} {self.id_type} {self.source}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,OrgID,IDType,Source,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType,Source]
        try:
            retInstance = cls.objects.get(organism_id=OrgID,id_type=IDType,source=Source)
        except:
            if verbose:
                print(f"[ID-Pub Not Found] {OrgID} {IDType} {Source}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgID,IDType,Source,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType,Source]
        return cls.objects.filter(organism_id=OrgID,id_type=IDType,source=Source).exists()

#=================================================================================================
class ID_Sequence(AuditModel):
    """
     Identification from Whole Genome Sequencing    
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "orgbatch_id":"OrgBatch ID",
        "id_type":"ID Type",
        "id_method":"ID Method",
        "id_organisms":"Organisms",
        "source": "Source",
        "id_date":"Date",
    }
    Choice_Dictionary = {
        'id_type':'ID_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    id_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "ID Method", on_delete=models.DO_NOTHING,
         db_column="id_type", related_name="%(class)s_idtype")
    run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID")
    id_method = models.CharField(max_length=25, blank=True, verbose_name = "Method")
    id_organisms =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Organisms", null=True, blank=True)
    id_date = models.DateField(blank=True, verbose_name = "ID Date")
    id_notes = models.CharField(max_length=120, blank=True,  verbose_name = "ID Notes")
    source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'id_seq'
        ordering=['orgbatch_id','id_type']
        indexes = [
             models.Index(name="idseq_drugid_idx",fields=['orgbatch_id']),
             models.Index(name="idseq_idtype_idx",fields=['id_type']),
             models.Index(name="idseq_idmed_idx",fields=['id_method']),
             models.Index(name="idseq_runid_idx",fields=['run_id']),
             models.Index(name="idseq_source_idx",fields=['source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.id_type} {self.run_id}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.id_type} {self.id_method} {self.run_id}"
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,IDType,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType,RunID]
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID,id_type=IDType,source=RunID)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {OrgBatchID} {IDType} {RunID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,IDType,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType,RunID]
        return cls.objects.filter(rgbatch_id=OrgBatchID,id_type=IDType,source=RunID).exists()

#=================================================================================================
class WGS_FastQC(AuditModel):
    """
     FastQC outcome - Import only   
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "orgbatch_id":"OrgBatch ID",
        "seq":"Seq",
        "run_id":"RunID",
        "base_stat" :"Statistics",
        "base_sequal" :"Per base sequence quality",
        "tile_sequal" :"Per tile sequence quality",
        "seq_qualsc" :"Per sequence quality scores",
        "base_gc" :"Per base sequence content",
        "seq_gc" :"Per sequence GC content",
        "base_N" :"Per base N content",
        "seq_len" :"Sequence Length Distribution",
        "seq_dup" :"Sequence Duplication Levels",
        "overrep" :"Overrepresented sequences",
        "adap_cont" :"Adapter Content",
    }
    Choice_Dictionary = {
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    seq = models.CharField(max_length=5, blank=True, verbose_name = "Seq")
    run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID") 
    base_stat = models.CharField(max_length=10, blank=True, verbose_name = "Basic Statistics")
    base_sequal = models.CharField(max_length=10, blank=True, verbose_name = "Per base sequence quality")
    tile_sequal = models.CharField(max_length=10, blank=True, verbose_name = "Per tile sequence quality")
    seq_qualsc = models.CharField(max_length=10, blank=True, verbose_name = "Per sequence quality scores")
    base_gc = models.CharField(max_length=10, blank=True, verbose_name = "Per base sequence content")
    seq_gc = models.CharField(max_length=10, blank=True, verbose_name = "Per sequence GC content")
    base_N = models.CharField(max_length=10, blank=True, verbose_name = "Per base N content")
    seq_len = models.CharField(max_length=10, blank=True, verbose_name = "Sequence Length Distribution")
    seq_dup = models.CharField(max_length=10, blank=True, verbose_name = "Sequence Duplication Levels")
    overrep = models.CharField(max_length=10, blank=True, verbose_name = "Overrepresented sequences")
    adap_cont = models.CharField(max_length=10, blank=True, verbose_name = "Adapter Content")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'wgs_fastqc'
        ordering=['orgbatch_id','seq','run_id']
        indexes = [
             models.Index(name="fastqc_orgbid_idx",fields=['orgbatch_id']),
             models.Index(name="fastqc_seq_idx",fields=['seq']),
             models.Index(name="fastqc_runid_idx",fields=['run_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.seq} {self.run_id}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.seq} {self.run_id}"
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,Seq,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID,Seq,RunID]
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID,seq=Seq,source=RunID)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {OrgBatchID} {Seq} {RunID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,Seq,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID,Seq,RunID]
        return cls.objects.filter(rgbatch_id=OrgBatchID,seq=Seq,source=RunID).exists()
    
#=================================================================================================
class WGS_CheckM(AuditModel):
    """
     CheckM outcome - Import only   
    """
#=================================================================================================
    HEADER_FIELDS   = {
        "orgbatch_id":"OrgBatch ID",
        "run_id":"RunID",
        "marker_lineage" :"Marker lineage",
        "n_genomes" :"n_genomes",
        "n_predit_genes" :"n_predit_genes",
        "n_markers" :"n_markers",
        "n_marker_sets" :"n_marker_sets",
        "n_scaffolds" :"n_scaffolds",
        "genome_size" :"genome_size",
        "coding_density" :"coding_density",
        "completeness" :"completeness",
        "contamination" :"contamination",
        "gc" : "gc",
        "gc_std" : "gc_std",
        "n_ambig_bases" :"n_ambig_bases",
        "longest_scaffold" :"longest_scaffold",
        "mean_scaffols" :"mean_scaffols",
        "n50" :"N50",
        "trans_table" :"trans_table",
    }
    Choice_Dictionary = {
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID")
    marker_lineage = models.CharField(max_length=25, blank=True, verbose_name = "Linage")
    n_genomes = models.IntegerField(default=0, blank=True, verbose_name ="n_genomes")
    n_predit_genes = models.IntegerField(default=0, blank=True, verbose_name ="n_predit_genes")
    n_markers = models.IntegerField(default=0, blank=True, verbose_name ="n_markers")
    n_marker_sets = models.IntegerField(default=0, blank=True, verbose_name ="n_marker_sets")
    n_scaffolds = models.IntegerField(default=0, blank=True, verbose_name ="n_scaffolds")
    genome_size = models.IntegerField(default=0, blank=True, verbose_name ="genome_size")
    coding_density = models.FloatField(default=0, blank=True, verbose_name ="coding_density")
    completeness = models.IntegerField(default=0, blank=True, verbose_name ="completeness")
    contamination = models.FloatField(default=0, blank=True, verbose_name ="contamination")
    gc = models.FloatField(default=0, blank=True, verbose_name ="gc")
    gc_std = models.FloatField(default=0, blank=True, verbose_name ="gc_std")
    n_ambig_bases = models.IntegerField(default=0, blank=True, verbose_name ="n_ambig_bases")
    longest_scaffold = models.IntegerField(default=0, blank=True, verbose_name ="longest_scaffold")
    mean_scaffols = models.FloatField(default=0, blank=True, verbose_name ="mean_scaffols")
    n50 = models.IntegerField(default=0, blank=True, verbose_name ="N50")
    trans_table = models.IntegerField(default=0, blank=True, verbose_name ="trans_table")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'wgs_checkm'
        ordering=['orgbatch_id','run_id']
        indexes = [
             models.Index(name="fastqc_orgbid_idx",fields=['orgbatch_id']),
             models.Index(name="fastqc_runid_idx",fields=['run_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.run_id}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.orgbatch_id} {self.run_id}"
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID,RunID]
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID,source=RunID)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {OrgBatchID} {Seq} {RunID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,RunID,verbose=0):
    # Returns an instance if found by [OrgBatchID,RunID]
        return cls.objects.filter(rgbatch_id=OrgBatchID,source=RunID).exists()