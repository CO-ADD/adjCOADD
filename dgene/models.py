from django.db import models

from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary
from dorganism.models import Taxonomy, Organism, Organism_Batch
from dscreen.models import Screen_Run


SEQUENCE_FILETYPES = [
        ('Reads','Raw reads'),
        ('FastQ','Fastq files (trimmed)'),
        ('FastA','Fasta files (assembly)'),
        ]

#=================================================================================================
# List of Sequences
#=================================================================================================
class Genome_Sequence(AuditModel):
#-------------------------------------------------------------------------------------------------
    LIST_VIEW_FIELDS = {
        #'seq_id':{"Seq ID":{"seq_id": URL_LINKS["seq_id"]},}, 
        'seq_id':"Seq ID", 
        "orgbatch_id.orgbatch_id":{'OrgBatch ID': {'orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "orgbatch_id.organism_id.organism_name":"Organism",
        'run_id':'Run ID',
        'seq_code':'SeqCode',  
        'seq_type':'Type',  
        'seq_method':'Method',  
        'runsample_file':'RunFile',
        'runsample_dir':'RunDir',
        'source':'Source',
    #    'source_code':'Source Code',
    #    'source_link':'Link',
    #    'reference':'Reference',
    #     'seq_date':'Seq Date'
        'seq_files':"Seq Files",
     }

    DICTIONARY_FIELDS = {
        'seq_type':'Seq_Type',      # DNA, RNA, cDNA, natDNA  
        'seq_method':'Seq_Method',  # Illumina, MinION-RB, MinION-Nat
    }

    ID_SEQUENCE = 'Sequence'
    ID_PREFIX = 'SEQ'
    ID_PAD = 5
    
    seq_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Seq ID")
    run_id = models.ForeignKey(Screen_Run, null=True, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 
    seq_code = models.CharField(max_length=120, unique=True, blank=True, verbose_name = "Seq Code")
    orgbatch_id = models.ForeignKey(Organism_Batch, null=True, blank=True, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id")

    seq_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Seq Type", on_delete=models.DO_NOTHING,
         db_column="seq_type", related_name="%(class)s_seqtype")
    seq_method = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Seq Method", on_delete=models.DO_NOTHING,
         db_column="seq_method", related_name="%(class)s_seqmethod")
    
    runsample_file = models.CharField(max_length=250, blank=True, verbose_name = "RunFile")
    runsample_dir = models.CharField(max_length=250, blank=True, verbose_name = "RunFolder")
    runsample_name = models.CharField(max_length=120, blank=True, verbose_name = "RunName")
     
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    source_link = models.CharField(max_length=120, blank=True, verbose_name = "Source Link")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")
#    seq_date = models.DateField(null=True, blank=True, verbose_name = "Seq Date")
    seq_files = models.CharField(max_length=50, blank=True, default='Reads', verbose_name = "Files")

    class Meta:
        app_label = 'dgene'
        db_table = 'genome_seq'
        ordering=['seq_id',]
        indexes = [
            models.Index(name="genoseq_seqid_idx",fields=['seq_id']),
            models.Index(name="genoseq_scr_idx",fields=['source']),
            models.Index(name="genoseq_seqcd_idx",fields=['seq_code']),
            models.Index(name="genoseq_scode_idx",fields=['source_code']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.seq_code}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.seq_code} ({self.seq_id}) "

   #------------------------------------------------
    @classmethod
    def get(cls,SeqID=None, SeqCode=None, verbose=0):
    # Returns an instance if found by [SeqID or SeqName]
        if SeqID:
            try:
                retInstance = cls.objects.get(seq_id=SeqID)
            except:
                if verbose:
                    print(f"[SeqID Not Found] {SeqID} ")
                retInstance = None
        elif SeqCode:
            try:
                retInstance = cls.objects.get(seq_code=SeqCode)
            except:
                if verbose:
                    print(f"[SeqCode Not Found] {SeqCode} ")
                retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,SeqID=None, SeqCode=None,verbose=0):
    # Returns an instance if found by [SeqID or SeqName]
        if SeqID:
            retValue = cls.objects.filter(seq_id=SeqID).exists()
        elif SeqCode:
            retValue = cls.objects.filter(seq_code=SeqCode).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.seq_id:
            self.seq_id = self.next_id()
            if self.seq_id: 
                super(Genome_Sequence, self).save(*args, **kwargs)
        else:
            super(Genome_Sequence, self).save(*args, **kwargs) 

    #--------------------------------------------------------------
    @staticmethod
    def gen_seq_code(OrgBatchID,RunID):
        return(f"{OrgBatchID}{SEQRUN_SEP}{RunID}")
    
#=================================================================================================
# Identification of Organism
#=================================================================================================
class ID_Pub(AuditModel):
    """
     Identification from Public or Collaborative sources    
    """
#=================================================================================================
    LIST_VIEW_FIELDS   = {
        "orgbatch_id.orgbatch_id":{'OrgBatch ID': {'orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "orgbatch_id.organism_id.organism_name":"Organism",
        "id_type":"ID Type",
        "id_method":"ID Method",
        "id_organisms":"Organisms",
        "source": "Source",
        "id_date":"Date",
        "id_notes":"Notes",
    }
    DICTIONARY_FIELDS = {
        'id_type':'ID_Type',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    # organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
    #     db_column="organism_id", related_name="%(class)s_organism_id") 
    id_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "ID Method", on_delete=models.DO_NOTHING,
         db_column="id_type", related_name="%(class)s_idtype")
    id_method = models.CharField(max_length=25, blank=True, verbose_name = "Method")
    id_organisms =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Organisms", null=True, blank=True)
    id_notes = models.CharField(max_length=120, blank=True,  verbose_name = "ID Notes")
    id_date = models.DateField(null=True, blank=True, verbose_name = "ID Date")
    source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'id_pub'
        ordering=['orgbatch_id','id_type']
        indexes = [
             models.Index(name="idp_org_idx",fields=['id_organisms']),
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
     Identification from Whole Genome Sequencing using Kraken, MLST and GTDBTK   
    """
#=================================================================================================
    LIST_VIEW_FIELDS   = {
        "seq_id.orgbatch_id.orgbatch_id":{'OrgBatch ID': {'seq_id.orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "seq_id.orgbatch_id.organism_id.organism_name":"Organism",
        "seq_id":"SeqID",
        "seq_id.run_id":'Run ID',
        "seq_filetype":"Seq Type", 
        "seq_file":"Seq File", 
        "kraken_organisms":"Kraken2 Organisms",
        "mlst_scheme": "MLST Scheme",
        "mlst_seqtype": "MLST SeqType",
        "mlst_alleles": "MLST Alleles",
        "gtdbtk_class": "GTDBTK",
        "gtdbtk_fastani": "FastANI",
        "source": "Source",
        "id_notes":"Notes",
   }
    DICTIONARY_FIELDS = {
        'seq_file':'Seq_File', # Trimmed, Contigs
        'seq_filetype':'Seq_FileType', # Trimmed, Contigs
    }


    seq_id = models.ForeignKey(Genome_Sequence, null=False, blank=False, verbose_name = "Seq ID", on_delete=models.DO_NOTHING,
        db_column="seq_id", related_name="%(class)s_seqid")
    seq_filetype =  models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Seq FType", on_delete=models.DO_NOTHING,
         db_column="seq_filetype", related_name="%(class)s_seqfiletype")
    seq_file = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Seq File", on_delete=models.DO_NOTHING,
         db_column="seq_file", related_name="%(class)s_seqfile")
    
    kraken_organisms =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Kraken2 Organisms", null=True, blank=True)
    mlst_scheme = models.CharField(max_length=20, blank=True, verbose_name = "MLST Scheme")
    mlst_seqtype = models.CharField(max_length=12, blank=True, verbose_name = "MLST SeqType")
    mlst_alleles = models.CharField(max_length=150, blank=True, verbose_name = "MLST Alleles")
    gtdbtk_class = models.CharField(max_length=120, blank=True, verbose_name = "MLST Scheme")
    gtdbtk_fastani = models.CharField(max_length=50, blank=True, verbose_name = "MLST SeqType")
    id_notes = models.CharField(max_length=120, blank=True,  verbose_name = "ID Notes")
    source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'id_seq'
        ordering=['seq_id__orgbatch_id','seq_id']
        indexes = [
#             models.Index(name="idseq_drugid_idx",fields=['orgbatch_id']),
             models.Index(name="idseq_seqfile_idx",fields=['seq_file']),
             models.Index(name="idseq_seqftype_idx",fields=['seq_filetype']),
             models.Index(name="idseq_seqid_idx",fields=['seq_id']),
             models.Index(name="idseq_source_idx",fields=['source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{str(self.seq_id)} {str(self.seq_filetype)} {self.seq_file} "
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{str(self.seq_id)} {str(self.seq_filetype)} {self.seq_file} "
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,SeqID,SeqFileType,SeqFile,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType,RunID]
        try:
            retInstance = cls.objects.get(seq_id=SeqID, seq_filetype=SeqFileType, seq_file=SeqFile,)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {SeqID} {SeqFile}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,SeqFile,SeqID,verbose=0):
    # Returns an instance if found by [OrgBatchID, IDType]
        return cls.objects.filter(seq_id__orgbatch_id=OrgBatchID,seq_file=SeqFile,seq_id=SeqID).exists()

#=================================================================================================
# Analysis of Whole Genome Sequence data
#=================================================================================================
class WGS_FastQC(AuditModel):
    """
     FastQC outcome from Trimming - Import only   
    """
#=================================================================================================
    LIST_VIEW_FIELDS   = {
        "seq_id.orgbatch_id.orgbatch_id":{'OrgBatch ID': {'seq_id.orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "seq_id.orgbatch_id.organism_id.organism_name":"Organism",
        "seq":"Seq",
        "seq_id":"SeqID",
        #"seq_id.run_id":'Run ID',
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
    DICTIONARY_FIELDS = {
    }

    seq = models.CharField(max_length=5, blank=True, verbose_name = "Seq")
    seq_id = models.ForeignKey(Genome_Sequence, null=False, blank=False, verbose_name = "Seq ID", on_delete=models.DO_NOTHING,
        db_column="seq_id", related_name="%(class)s_seqid") 
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
        ordering=['seq_id__orgbatch_id','seq','seq_id']
        indexes = [
        #     models.Index(name="fastqc_orgbid_idx",fields=['orgbatch_id']),
             models.Index(name="fastqc_seq_idx",fields=['seq']),
             models.Index(name="fastqc_seqid_idx",fields=['seq_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.seq_id__orgbatch_id} {self.seq} {str(self.seq_id)}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.seq_id__orgbatch_id} {self.seq} {str(self.seq_id)}"
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,Seq,SeqID,verbose=0):
    # Returns an instance if found by [OrgBatchID,Seq,RunID]
        try:
            retInstance = cls.objects.get(seq_id__orgbatch_id=OrgBatchID,seq=Seq,seq_id=SeqID)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {OrgBatchID} {Seq} {SeqID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,Seq,SeqID,verbose=0):
    # Returns an instance if found by [OrgBatchID,Seq,RunID]
        return cls.objects.filter(seq_id__orgbatch_id=OrgBatchID,seq=Seq,seq_id=SeqID).exists()
    
#=================================================================================================
class WGS_CheckM(AuditModel):
    """
     CheckM outcome of Assemblies- Import only   
    """
#=================================================================================================
    LIST_VIEW_FIELDS   = {
        "seq_id.orgbatch_id.orgbatch_id":{'OrgBatch ID': {'seq_id.orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "seq_id.orgbatch_id.organism_id.organism_name":"Organism",
        "seq_id":"SeqID",
        #"seq_id.run_id":'Run ID',
        "assembly":"Assembly",
        "fasta":"FastA",
        "assembly_qc":"QC",
        "marker_lineage" :"Marker lineage",
        "completeness" :"Completeness",
        "contamination" :"Contamination",
        "n_genomes" :"#Genomes",
        "n_predit_genes" :"#Predit Genes",
        "n_markers" :"#Markers",
        "n_marker_sets" :"#Marker Sets",
        "genome_size" :"Genome Size",
        "coding_density" :"Coding Density",
        "n_contigs" :"#Contigs",
        "longest_contig" :"Longest Contig",
        "mean_contigs" :"Mean Contigs",
        "n50_contigs" :"N50 Contigs",
        "gc" : "GC",
        "gc_std" : "GC Std",
        "n_ambig_bases" :"#Ambig Bases",
        "trans_table" :"Trans Table",
    }
    DICTIONARY_FIELDS = {
    }

    seq_id = models.ForeignKey(Genome_Sequence, null=False, blank=False, verbose_name = "Seq ID", on_delete=models.DO_NOTHING,
        db_column="seq_id", related_name="%(class)s_seqid")
    assembly = models.CharField(max_length=25, blank=True, verbose_name = "Assembly")
    fasta = models.CharField(max_length=55, blank=True, verbose_name = "Fasta")
    assembly_qc = models.CharField(max_length=15, blank=True, verbose_name = "Assembly QC") 
    marker_lineage = models.CharField(max_length=25, blank=True, verbose_name = "Linage")
    n_genomes = models.IntegerField(default=0, blank=True, verbose_name ="n_genomes")
    n_predit_genes = models.IntegerField(default=0, blank=True, verbose_name ="n_predit_genes")
    n_markers = models.IntegerField(default=0, blank=True, verbose_name ="n_markers")
    n_marker_sets = models.IntegerField(default=0, blank=True, verbose_name ="n_marker_sets")
    n_contigs = models.IntegerField(default=0, blank=True, verbose_name ="n_contigs")
    genome_size = models.IntegerField(default=0, blank=True, verbose_name ="genome_size")
    
    coding_density = models.DecimalField(max_digits=9, decimal_places=3, default=0, blank=True, verbose_name ="coding_density")
    completeness = models.DecimalField(max_digits=9, decimal_places=2, default=0, blank=True, verbose_name ="completeness")
    contamination = models.DecimalField(max_digits=9, decimal_places=2, default=0, blank=True, verbose_name ="contamination")
    gc = models.DecimalField(max_digits=9, decimal_places=3,default=0, blank=True, verbose_name ="gc")
    gc_std = models.DecimalField(max_digits=9, decimal_places=3,default=0, blank=True, verbose_name ="gc_std")
    n_ambig_bases = models.IntegerField(default=0, blank=True, verbose_name ="n_ambig_bases")
    longest_contig = models.IntegerField(default=0, blank=True, verbose_name ="longest_contig")
    mean_contigs = models.DecimalField(max_digits=12, decimal_places=2, default=0, blank=True, verbose_name ="mean_contigs")
    n50_contigs = models.IntegerField(default=0, blank=True, verbose_name ="N50_contigs")
    trans_table = models.IntegerField(default=0, blank=True, verbose_name ="trans_table")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'wgs_checkm'
        ordering=['seq_id__orgbatch_id','seq_id']
        indexes = [
        #     models.Index(name="checkqc_orgbid_idx",fields=['orgbatch_id']),
             models.Index(name="checkqc_seqid_idx",fields=['seq_id']),
             models.Index(name="checkqc_ass_idx",fields=['assembly']),
             models.Index(name="checkqc_fa_idx",fields=['fasta']),
             models.Index(name="checkqc_comp_idx",fields=['completeness']),
             models.Index(name="checkqc_cont_idx",fields=['contamination']),
             models.Index(name="checkqc_qc_idx",fields=['assembly_qc']),
             models.Index(name="checkqc_ml_idx",fields=['marker_lineage']),
        ]

    #------------------------------------------------
    # def __str__(self) -> str:
    #     retStr = f"{self.orgbatch_id} {str(self.seq_id)}"
    #     return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.seq_id} [{self.assembly} - {self.fasta}]"
        return(retStr)


   #------------------------------------------------
    @classmethod
    def get(cls,SeqID,Assembly,FastA,verbose=0):
    # Returns an instance if found by [OrgBatchID,RunID]
        try:
            retInstance = cls.objects.get(seq_id=SeqID,assembly=Assembly,fasta=FastA)
        except:
            if verbose:
                print(f"[ID-WGS Not Found] {SeqID} - {Assembly} - {FastA}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,SeqID,Assembly,verbose=0):
    # Returns an instance if found by [OrgBatchID,RunID]
        return cls.objects.filter(orgbatch_id=OrgBatchID,seq_id=SeqID,assembly=Assembly).exists()
    
#=================================================================================================
# Gene of Interest
#=================================================================================================
class Gene(AuditModel):
    """
    List of Genes
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
        #"gene_id":{"Gene Name":{"gene_id": URL_LINKS["gene_id"]},},
        "gene_id":"Gene ID",
        "gene_code":"Gene Code",
        "gene_note":"Gene Note",
        "gene_type":"Gene Type",
        "gene_subtype":"Gene SubType",
        "amr_class":"AMR Class",
        "amr_subclass":"AMR SubClass",
        "protein_class":"Protein Class",
        "source": "Source",
    }

    DICTIONARY_FIELDS = {
        'gene_type':'Gene_Type',
    }

    FORM_GROUPS={
       'Group_gene': ["source", "gene_note", "gene_type", "organisms", "resistance_class", "uniprot", "variants"],
       'Group_protein': ['protein_name','protein_class','protein_subclass','gene_note','gene_modification'],
    }

    ID_SEQUENCE = 'Gene'
    ID_PREFIX = 'GEN'
    ID_PAD = 5

    gene_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Gene ID")
    gene_code = models.CharField(max_length=25, unique=True,  verbose_name = "Gene Code")
    gene_name = models.CharField(max_length=25, blank=True,  verbose_name = "Gene Name")
    gene_note = models.CharField(max_length=500, verbose_name = "Gene Note")
    urlname = models.SlugField(max_length=30, verbose_name = "URLGene")
    gene_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Gene Type", on_delete=models.DO_NOTHING,
         db_column="gene_type", related_name="%(class)s_genetype")
    gene_subtype = models.CharField(max_length=50, blank=True,  verbose_name = "Gene Subtype")
    amr_class = models.CharField(max_length=50, blank=True,  verbose_name = "AMR Class")
    amr_subclass = models.CharField(max_length=250, blank=True,  verbose_name = "AMR Subclass")
    gene_modification = models.CharField(max_length=120, verbose_name = "Modification")
    protein_name = models.CharField(max_length=150, blank=True,  verbose_name = "Protein")
    protein_class = models.CharField(max_length=50, blank=True,  verbose_name = "Class")
    protein_subclass = models.CharField(max_length=50,  blank=True, verbose_name = "SubClass")
    resistance_class = models.CharField(max_length=50, blank=True,  verbose_name = "Resistance Class")
    variants = models.CharField(max_length=120, blank=True,  verbose_name = "Variants")
    genbank = models.CharField(max_length=120, blank=True,  verbose_name = "GenBank")
    uniprot = models.CharField(max_length=120, blank=True,  verbose_name = "UniProt")
    organisms =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=10, verbose_name = "Organisms", null=True, blank=True)
    source = models.CharField(max_length=150, blank=True, verbose_name = "Source")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'gene'
        ordering=['gene_code','amr_class']
        indexes = [
             models.Index(name="gene_gtype_idx",fields=['gene_type']),
             models.Index(name="gene_aclass_idx",fields=['amr_class']),
             models.Index(name="gene_asclass_idx",fields=['amr_subclass']),
             models.Index(name="gene_gmod_idx",fields=['gene_modification']),
             models.Index(name="gene_pclass_idx",fields=['protein_class']),
             models.Index(name="gene_psclass_idx",fields=['protein_subclass']),
             models.Index(name="gene_src_idx",fields=['source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.gene_code}"
        return(retStr)

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.gene_code} {self.gene_type}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,GeneID=None, GeneCode=None, verbose=0):
    # Returns an instance if found by [SeGeneIDqID or GeneCode]
        if GeneID:
            try:
                retInstance = cls.objects.get(gene_id=GeneID)
            except:
                if verbose:
                    print(f"[GeneID Not Found] {GeneID} ")
                retInstance = None
        elif GeneCode:
            try:
                retInstance = cls.objects.get(gene_code=GeneCode)
            except:
                if verbose:
                    print(f"[GeneCode Not Found] {GeneCode} ")
                retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,GeneID=None, GeneCode=None,verbose=0):
    # Returns an instance if found by [GeneID or GeneCode]
        if GeneID:
            retValue = cls.objects.filter(gene_id=GeneID).exists()
        elif GeneCode:
            retValue = cls.objects.filter(gene_code=GeneCode).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        self.urlname = slugify(self.gene_code,allow_unicode=False)
        if not self.gene_id:
            self.gene_id = self.next_id()
            if self.gene_id: 
                super(Gene, self).save(*args, **kwargs)
        else:
            super(Gene, self).save(*args, **kwargs) 

#=================================================================================================
class AMR_Genotype(AuditModel):
    """
    List of AMR Genes
    """
#=================================================================================================

    LIST_VIEW_FIELDS = {
        "seq_id.orgbatch_id.orgbatch_id":{'OrgBatch ID': {'seq_id.orgbatch_id.organism_id.organism_id':URL_LINKS["organism_id"]}},
        "seq_id.orgbatch_id.organism_id.organism_name":"Organism",
        "seq_id":"SeqID",
        "seq_id.run_id":'Run ID',
        #"gene_id":{"Gene Name":{"gene_id": URL_LINKS["gene_id"]},},
        "gene_id.gene_code":"Gene Code",
        "gene_id.gene_type":"Gene Type",
        "gene_id.amr_class":"AMR Class",
        "gene_id.amr_subclass":"AMR SubClass",
        #"amr_method":"Method",
        #"seq_method":"Seq Method",
        "seq_coverage":"Coverage",
        "seq_identity":"Identity",
        "closest_id": "Closest",
        "closest_name": "Closest Name",
    }

    seq_id = models.ForeignKey(Genome_Sequence, null=True, blank=True, verbose_name = "Seq ID", on_delete=models.DO_NOTHING,
        db_column="seq_id", related_name="%(class)s_seqid")
    gene_id = models.ForeignKey(Gene, null=False, blank=False, verbose_name = "Gene ID", on_delete=models.DO_NOTHING,
        db_column="gene_id", related_name="%(class)s_geneid")

    amr_method = models.CharField(max_length=25, blank=True,   verbose_name = "AMR Method")

    seq_method = models.CharField(max_length=25, blank=True,    verbose_name = "Seq Method")
    seq_coverage = models.DecimalField(max_digits=9, decimal_places=2,default=0, blank=True, verbose_name ="Seq Coverage") 
    seq_identity = models.DecimalField(max_digits=9, decimal_places=2,default=0, blank=True, verbose_name ="Seq Identity")
    closest_id = models.CharField(max_length=150, blank=True,    verbose_name = "Closest")
    closest_name = models.CharField(max_length=150, blank=True, verbose_name = "Closest Name")

    #------------------------------------------------
    class Meta:
        app_label = 'dgene'
        db_table = 'amr_genotype'
        ordering=['gene_id','amr_method','seq_id']
        indexes = [
             models.Index(name="amrgt_am_idx",fields=['amr_method']),
             models.Index(name="amrgt_gid_idx",fields=['gene_id']),
             models.Index(name="amrgt_sid_idx",fields=['seq_id']),
        #     models.Index(name="amrgt_obid_idx",fields=['orgbatch_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        retStr = f"{self.pk} {self.gene_id.gene_code} {self.seq_id.orgbatch_id} {self.seq_id}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,GeneID, Method, SeqID=None, OrgBatchID=None, verbose=0):
    # Returns an instance if found by [SeGeneIDqID or GeneCode]
        if SeqID:
            try:
                retInstance = cls.objects.get(gene_id=GeneID,seq_id=SeqID,amr_method=Method)
            except:
                if verbose:
                    print(f"[AMR_Genotype Not Found] {GeneID} {SeqID} {Method} ")
                retInstance = None
        elif OrgBatchID:
            try:
                retInstance = cls.objects.get(gene_id=GeneID,seq_id__orgbatch_id=OrgBatchID,amr_method=Method)
            except:
                if verbose:
                    print(f"[AMR_Genotype Not Found] {GeneID} {OrgBatchID} {Method} ")
                retInstance = None
        else:
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,GeneID, Method, SeqID=None, OrgBatchID=None,verbose=0):
    # Returns an instance if found by [GeneID or GeneCode]
        if GeneID:
            retValue = cls.objects.filter(gene_id=GeneID,seq_id=SeqID,amr_method=Method).exists()
        elif OrgBatchID:
            retValue = cls.objects.filter(gene_id=GeneID,seq_id__orgbatch_id=OrgBatchID,amr_method=Method).exists()
        else:
            retValue = False
        return(retValue)
