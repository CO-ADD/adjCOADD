import re

#from django_rdkit import models
from django.db import models
from django_rdkit import models as rdm
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex, GistIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from applib.bio.biln import sync_sequences, fasta_to_biln, biln_residue_count, biln_molecular_weight
from dorganism.models import Taxonomy 
from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

# Natural amino-acid one-letter codes (FASTA).
_NATURAL_AA = set('ACDEFGHIKLMNPQRSTVWY')

#-- TO DO -------------------------------------------------------------
"""
    Remove Biologic from Isoform 
    Modification_Dictionary
        Add Description, Cleavage sites ?
    Isoform_Batch
        Add volume, concentration, purity, storage conditions
        
    Sequence properties
        MW, nSeq, pI, hydrophobicity, charge, 
"""
#-----------------------------------------------------------------------

    
def _sync_sequence_fields(instance):
    """Detect edited sequence fields and synchronise biln/helm/fasta.

    Only conversions triggered by an edited field run. For a new instance all
    three fields are treated as new. BILN is the master source when several
    fields change at once (priority BILN > HELM > FASTA).
    """
    changed = set()
    if instance.pk:
        db = type(instance).objects.filter(pk=instance.pk).values_list(
            'biln', 'helm', 'fasta').first()
        if db is not None:
            db_biln, db_helm, db_fasta = (v or '' for v in db)
            if (instance.biln or '') != db_biln:
                changed.add('biln')
            if (instance.helm or '') != db_helm:
                changed.add('helm')
            if (instance.fasta or '') != db_fasta:
                changed.add('fasta')
        else:
            changed = {'biln', 'helm', 'fasta'}
    else:
        changed = {'biln', 'helm', 'fasta'}

    instance.biln, instance.helm, instance.fasta = sync_sequences(
        instance.biln, instance.helm, instance.fasta, changed)


def _compute_mw_nseq(instance):
    """Update molecular weight (mw) and residue count (nseq) from the sequence.

    nseq is the amino-acid count from FASTA when available, falling back to the
    number of residues in the main BILN chain; mw is the average molecular
    weight of the RDKit molecule built from BILN. All failures degrade
    gracefully (values left at 0).
    """
    biln = (instance.biln or '').strip()
    fasta = (instance.fasta or '').strip()

    if fasta:
        instance.nseq = len(fasta)
    elif biln:
        try:
            instance.nseq = biln_residue_count(biln)
        except Exception:
            instance.nseq = 0
    else:
        instance.nseq = 0

    if biln:
        try:
            mw = biln_molecular_weight(biln)
            instance.mw = round(mw, 3) if mw else 0
        except Exception:
            instance.mw = 0
    else:
        instance.mw = 0


def _compute_sequence_properties(instance):
    """Update sequence-derived properties (mw, nseq, pI, gravy, charge).

    MW/nseq via ``_compute_mw_nseq``; isoelectric_point / gravy / charge_at_ph7
    are computed from the FASTA sequence via BioPython and left NULL for
    modified/non-natural sequences. All failures degrade gracefully.
    """
    _compute_mw_nseq(instance)

    fasta = (instance.fasta or '').strip()

    # pI / hydrophobicity / charge (BioPython, natural sequences only)
    pa = _protein_analysis(fasta)
    if pa is not None:
        try:
            instance.isoelectric_point = round(pa.isoelectric_point(), 2)
        except Exception:
            instance.isoelectric_point = None
        try:
            instance.gravy = round(pa.gravy(), 3)
        except Exception:
            instance.gravy = None
        try:
            instance.charge_at_ph7 = round(pa.charge_at_pH(7.0), 2)
        except Exception:
            instance.charge_at_ph7 = None
    else:
        instance.isoelectric_point = None
        instance.gravy = None
        instance.charge_at_ph7 = None


def _protein_analysis(fasta):
    """Return a BioPython ProteinAnalysis for a natural FASTA, or None.

    Modified/non-natural sequences (and any BioPython import failure) return
    None so callers can fall back gracefully.
    """
    fasta = ''.join(str(fasta or '').split()).upper()
    if not fasta or any(ch not in _NATURAL_AA for ch in fasta):
        return None
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        return ProteinAnalysis(fasta)
    except Exception:
        return None


def _sequence_segment(instance):
    """Return the BILN segment for a Biologic / Modification_Dictionary instance.

    Prefers the stored BILN; otherwise derives BILN from FASTA. Empty on
    failure so concatenation never inserts a broken segment.
    """
    biln = (getattr(instance, 'biln', '') or '').strip()
    if biln:
        return biln
    fasta = (getattr(instance, 'fasta', '') or '').strip()
    if fasta:
        try:
            return fasta_to_biln(fasta)
        except Exception:
            return ''
    return ''


#-------------------------------------------------------------------------------------------------
# Biologic Application Model
# Biologic
#  -> Biologic_Isoform  <- Modification_Dictionary (via Biologic_Modification)
#     -> Isoform_Batch
#
# AminoAcid (dchem.models) - List of Amino Acids for Biologic Sequences, including modified and non-natural amino acids
# Modification_Dictionary - List of Biologic Modifications/Tags, including modified and non-natural amino acids
#

#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Biologic(AuditModel):
    """
    Main class of Biologics s
    
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
#       'organism_name':{"VerboseName":'Biologic Name','Updatable':False}
        'biologic_id':{'Biologic ID': {'biologic_id':URL_LINKS['biologic_id']}}, 
        'biologic_name':'Name',
        'biologic_class':'Type',
        # 'biologic_panel':'Panel',
        'biologic_notes':'Notes',
        # 'biologic_code':'Biologic Code',
        # 'biologic_identification':'Identification',
        # 'biologic_origin':'Origin',
        'biln':'Biologic Sequence',
        # 'source':"Source",
        # 'source_code':"Source Code",
        # 'reference': "Reference",
        #'tax_id':{'Tax-ID': {'tax_id':URL_LINKS['tax_id']}}, 
    }

    CARDS_FIELDS= {
        "source" : "Source",
        "biologic_class": "Class",
    #     "res_property": "Phenotype",
    #     "gen_property": "Genotype",
    }

    DETAIL_VIEW_GROUPS = [
        ['biologic_name', 'biologic_class', 'biologic_panel', 'biologic_notes', 'biologic_origin','source', 'source_code','reference','biologist'],
        ['mw', 'nseq'],
        ['biln',],
        #['mta_status','mta_notes','mta_document','biologist'],
    ]

    DICTIONARY_FIELDS = {
        'mta_status':'License_Status',
        'biologic_class':'Peptide_Class',
        'biologic_panel':'Peptide_Panel',
    }

    biologic_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Biologic ID") 
    biologic_name= models.CharField(max_length=200, blank=True, verbose_name = "Biologic Name") 
    biologic_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Biologic Notes")
    #biologic_code= models.CharField(max_length=30, blank=True, verbose_name = "Biologic Code")
    biologic_panel=ArrayField(models.CharField(max_length=100, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    biologic_class=models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="biologic_class", related_name="%(class)s_biologic_class")
    #biologic_identification = models.CharField(max_length=512, blank=True, verbose_name = "Biologic Identification")
    biologic_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Biologic")

    biln = models.TextField(blank=True, verbose_name = "Sequence (BILN)") 
    helm = models.TextField(blank=True, verbose_name = "Sequence (HELM)") 
    fasta = models.TextField(blank=True, verbose_name = "Sequence (FASTA)") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    nseq = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nSeq")
    isoelectric_point = models.DecimalField(max_digits=6, decimal_places=2, null=True, blank=True, verbose_name ="pI")
    gravy = models.DecimalField(max_digits=8, decimal_places=3, null=True, blank=True, verbose_name ="Hydrophobicity (GRAVY)")
    charge_at_ph7 = models.DecimalField(max_digits=6, decimal_places=2, null=True, blank=True, verbose_name ="Charge @ pH 7")

    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    # organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
    #     db_column="organism_name", related_name="%(class)s_organism_name")

    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_mta")
    mta_notes = models.CharField(max_length=512, blank=True, verbose_name='MTA Notes')
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")

    # collect_tissue = models.CharField(max_length=120, blank=True, verbose_name = "From Tissue/Organ")
    # patient_diagnosis = models.CharField(max_length=120, blank=True, verbose_name = "Patient Diagnosis")
    # patient = models.CharField(max_length=20, blank=True, verbose_name = "Patient Info")

    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    assoc_documents = models.ManyToManyField(Document,verbose_name = "Documents", blank=True,
        db_table = "biologic_doc", related_name="%(class)s_document")
    
    #tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")


#------------------------------------------------

    class Meta:
        app_label = 'dbiologic'
        db_table = 'biologic'
        #ordering=['biologic_id']
        indexes = [
            #models.Index(name="pep_code_idx", fields=['biologic_code']),
            models.Index(name="pep_type_idx", fields=['biologic_class']),
            GinIndex(name="pep_panel_idx", fields=['biologic_panel']),
            #models.Index(name="pep_source_idx", fields=['source']),
            # models.Index(name="org_taxid_idx", fields=['tax_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.biologic_id} ({self.biologic_name})"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.biologic_id} {self.biologic_name}".strip()

    #------------------------------------------------
    @classmethod
    def exists(cls,BiologicID=None,BiologicName=None,verbose=0):
        if BiologicID:
            # Returns if an instance exists by biologic_id
            return cls.objects.filter(biologic_id=BiologicID).exists()
        elif BiologicName:
            # Returns if an instance exists by biologic_line
            return cls.objects.filter(biologic_name=BiologicName).exists()

    #------------------------------------------------
    @classmethod
    def get(cls,BiologicID=None,BiologicName=None,verbose=0):
        if BiologicID:
            # Returns an instance by biologic_id
            try:
                retInstance = cls.objects.get(biologic_id=BiologicID)
            except:
                if verbose:
                    print(f"[BiologicID Not Found] {BiologicID} ")
                retInstance = None
        elif BiologicName:
            # Returns an instance by biologic_id
            try:
                retInstance = cls.objects.get(biologic_name=BiologicName)
            except:
                if verbose:
                    print(f"[Biologic Name Not Found] {BiologicName} ")
                retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def str_BiologicID(cls,BiologicClass,BiologicNo) -> str:
    #
    # Output:   Biologic_ID as string like MAB_0001, PEP_0001, NB_0001 
    #
        return(f"{BiologicClass}{PEPID_SEP}{BiologicNo:04d}")


    #------------------------------------------------
    @classmethod
    def find_Next_BiologicID(cls,BiologicClass, BiologicClassType = PEPTIDE_CLASSES) -> str:
        if BiologicClass in BiologicClassType:
            Biologic_IDSq=Sequence(BiologicClass)
            Biologic_nextID = next(Biologic_IDSq)
            Biologic_strID = cls.str_BiologicID(BiologicClass,Biologic_nextID)
            while cls.exists(Biologic_strID):
                Biologic_nextID = next(Biologic_IDSq)
                Biologic_strID = cls.str_BiologicID(BiologicClass,Biologic_nextID)
            return(Biologic_strID)    
        else:
            return(None)

     #------------------------------------------------
    def save(self, *args, **kwargs):
        _sync_sequence_fields(self)
        _compute_sequence_properties(self)
        if not self.biologic_id:
            _class = str(self.biologic_class.dict_value) if self.biologic_class else 'PEP'
            self.biologic_id = self.find_Next_BiologicID(_class)
            if self.biologic_id:
                super(Biologic, self).save(*args, **kwargs)
        else:
            super(Biologic, self).save(*args, **kwargs) 


#=================================================================================================
class Biologic_Isoform(AuditModel):
    """
    Isoform - a modified variant of a Biologic
    
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
        'isoform_id': {'Isoform ID': {'isoform_id': URL_LINKS['isoform_id']}},
        'isoform_name': 'Isoform Name',
        'biologic_id': {'Biologic ID': {'biologic_id': URL_LINKS['biologic_id']}},
        #'biln': 'Sequence',
        'mw': 'MW',
        'biologist': 'Biologist',
    }

    DETAIL_VIEW_GROUPS = [
        ['isoform_id', 'isoform_name', 'isoform_notes', 'biologic_id', 'mw', 'nseq', 'biologist'],
        ['biln'],
        #['biologic_modifications'],
    ]
    DICTIONARY_FIELDS = {}

    isoform_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Isoform ID") 
    isoform_name= models.CharField(max_length=200, blank=True, verbose_name = "Isoform Name") 
    isoform_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Isoform Notes")
    
    biologic_id = models.ForeignKey(Biologic, null=False, blank=False, verbose_name = "Biologic ID", on_delete=models.DO_NOTHING,
        db_column="biologic_id", related_name="isoforms")

    biologic_modifications = models.ManyToManyField(
        'Modification_Dictionary',
        through='Biologic_Modification',
        through_fields=('isoform_id', 'biomod_id'),
        blank=True,
        related_name='isoforms',
        verbose_name = "Modifications",
    )

    biln = models.TextField(blank=True, verbose_name = "Sequence (BILN)") 
    helm = models.TextField(blank=True, verbose_name = "Sequence (HELM)") 
    fasta = models.TextField(blank=True, verbose_name = "Sequence (FASTA)") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    nseq = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nSeq")
    isoelectric_point = models.DecimalField(max_digits=6, decimal_places=2, null=True, blank=True, verbose_name ="pI")
    gravy = models.DecimalField(max_digits=8, decimal_places=3, null=True, blank=True, verbose_name ="Hydrophobicity (GRAVY)")
    charge_at_ph7 = models.DecimalField(max_digits=6, decimal_places=2, null=True, blank=True, verbose_name ="Charge @ pH 7")

    assoc_documents = models.ManyToManyField(Document,verbose_name = "Documents", blank=True,
        db_table = "isoform_doc", related_name="%(class)s_document")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    class Meta:
        app_label = 'dbiologic'
        db_table = 'biologic_isoform'
        ordering=['biologic_id','isoform_id']
        indexes = [
            models.Index(name="pepiso_pep_idx", fields=['biologic_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.isoform_id} ({self.isoform_name})"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.isoform_id} {self.isoform_name}".strip()

    #------------------------------------------------
    @classmethod
    def str_IsoformID(cls, BiologicID, IsoformNo) -> str:
        return(f"{BiologicID}{PEPISO_SEP}I{IsoformNo:02d}")

    #------------------------------------------------
    def find_Next_IsoformID(self) -> str:
        biologic_id = self.biologic_id.biologic_id
        next_no = 1
        next_iso = self.str_IsoformID(biologic_id, next_no)
        while self.exists(next_iso):
            next_no += 1
            next_iso = self.str_IsoformID(biologic_id, next_no)
        return(next_iso)

    #------------------------------------------------
    @classmethod
    def get(cls, IsoformID, verbose=0):
        try:
            retInstance = cls.objects.get(isoform_id=IsoformID)
        except:
            if verbose:
                print(f"[Isoform Not Found] {IsoformID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls, IsoformID, verbose=0):
        return cls.objects.filter(isoform_id=IsoformID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        _sync_sequence_fields(self)
        _compute_sequence_properties(self)
        if not self.isoform_id:
            self.isoform_id = self.find_Next_IsoformID()
            if self.isoform_id:
                super(Biologic_Isoform, self).save(*args, **kwargs)
        else:
            super(Biologic_Isoform, self).save(*args, **kwargs)

    #------------------------------------------------
    def rebuild_sequence_from_modifications(self, *args, **kwargs):
        """Rebuild the isoform sequence from the parent biologic plus its
        terminal modifications, then re-sync and re-save.

        Assembly order (per project convention):
            C-terminal modifications (ordered by mod_order) at the FRONT,
            then the parent biologic sequence,
            then N-terminal modifications (ordered by mod_order) at the END.
        Internal ('I') modifications are not concatenated.
        """
        parts = []
        for link in self.modification_links.filter(mod_position='C').order_by('mod_order'):
            seg = _sequence_segment(link.biomod_id)
            if seg:
                parts.append(seg)

        if self.biologic_id:
            seg = _sequence_segment(self.biologic_id)
            if seg:
                parts.append(seg)

        for link in self.modification_links.filter(mod_position='N').order_by('mod_order'):
            seg = _sequence_segment(link.biomod_id)
            if seg:
                parts.append(seg)

        self.biln = '-'.join(parts)
        self.save(*args, **kwargs)

    #------------------------------------------------
    @property
    def nterm_modifications(self):
        return self.modification_links.filter(mod_position='N').order_by('mod_order')

    @property
    def cterm_modifications(self):
        return self.modification_links.filter(mod_position='C').order_by('mod_order')

    @property
    def internal_modifications(self):
        return self.modification_links.filter(mod_position='I').order_by('mod_order')

#=================================================================================================
class Isoform_Batch(AuditModel):
    """
    Biologic Isoform Batch Collection
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
        "batch_id":"Batch ID",
        "batch_notes":"Batch Notes",
        "source":"Source",
        "source_code":"Code",
        "source_type":"Method",
        "qc_status":"QC",
        "quality_source": "Quality by",
        "conc":"Conc",
        "conc_unit":"Conc Units",
        "init_amount":"Amount",
        "init_amount_unit":"Amount Unit",
        "biologist":"Biologist"
    }

    DICTIONARY_FIELDS = {
        'qc_status':'QC_Status',
        'source_type' : 'Peptide_Source',
        'conc_unit' : 'Concentration_Units',
        'init_amount_unit' : 'Amount_Units',
    }
    
    FORM_GROUPS = {
       'Group1': ["batch_id", "batch_notes", "qc_status", "quality_source", "stock_date", "stock_level", "biologist", "source", "source_code", "source_type", "conc", "conc_unit", "init_amount", "init_amount_unit" ]
       }

    isoformbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "IsoFormBatch ID")
    isoform_id = models.ForeignKey(Biologic_Isoform, null=False, blank=False, verbose_name = "Isoform ID", on_delete=models.DO_NOTHING,
        db_column="isoform_id", related_name="batches")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    source_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Source Type", on_delete=models.DO_NOTHING,
        db_column="source_type", related_name="%(class)s_source_type")
    
    conc = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="Concentration")
    conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Concentration Unit", on_delete=models.DO_NOTHING,
        db_column="conc_unit", related_name="%(class)s_conc_unit")
    
    init_amount = models.DecimalField(max_digits=12, decimal_places=3, default=0, verbose_name = "Init Amount")
    init_amount_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Init Amount Unit", on_delete=models.DO_NOTHING,
        db_column="init_amount_unit", related_name="%(class)s_init_amount_unit")

    quality_source = models.CharField(max_length=150, blank=True, verbose_name = "QC Source")
    qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC status", on_delete=models.DO_NOTHING,
        db_column="qc_status", related_name="%(class)s_qc")
    qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")

    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date") 
    stock_level = models.CharField(max_length=20, blank=True, verbose_name = "Stock Levels") 
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dbiologic'
        db_table = 'isoform_batch'
        ordering=['isoformbatch_id']
        indexes = [
            models.Index(name="isobatch_isobatch_idx",fields=['isoform_id','batch_id']),
            models.Index(name="isobatch_qc_idx",fields=['qc_status']),
            models.Index(name="isobatch_sdate_idx",fields=['stock_date']),
            models.Index(name="isobatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.isoformbatch_id}"

    #------------------------------------------------
    @classmethod
    # Formats BatchNo:int -> BatchID:str 
    def str_BatchID(self,BatchNo:int) -> str:
        return(f"{BatchNo:02d}")
    #------------------------------------------------
    @classmethod
    # Formats IsoformID:str,BatchID:str -> IsoBatchID:str
    def str_IsoBatchID(self,IsoformID:str,BatchID:str) -> str:
        return(f"{IsoformID}{PEPBATCH_SEP}{BatchID}")

    #------------------------------------------------
    def find_Next_BatchID(self, IsoformID:str, BatchID:str=None) -> str:
        # Check for given BatchID    
        if BatchID:
            # Clean up BatchID - remove non alphanumeric character and make uppercase
            BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

            # Clean up BatchID - reformat numbers
            if BatchID.isnumeric():
                BatchID = self.str_BatchID(int(BatchID))

            if not self.exists(self.str_IsoBatchID(IsoformID, BatchID)):
                return(BatchID)

        # Find new BatchID    
        next_BatchNo = 1
        next_IsoBatch = self.str_IsoBatchID(IsoformID, self.str_BatchID(next_BatchNo))
        while self.exists(next_IsoBatch):
            next_BatchNo = next_BatchNo + 1
            next_IsoBatch = self.str_IsoBatchID(IsoformID, self.str_BatchID(next_BatchNo))
        return(self.str_BatchID(next_BatchNo))    

    #------------------------------------------------
    @classmethod
    def get(cls,IsoBatchID,verbose=0):
    # Returns an instance if found by isoformbatch_id
        try:
            retInstance = cls.objects.get(isoformbatch_id=IsoBatchID)
        except:
            if verbose:
                print(f"[IsoBatch Not Found] {IsoBatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,IsoBatchID,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(isoformbatch_id=IsoBatchID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.isoformbatch_id: 
            # creates new IsoBatchID
            IsoformID = self.isoform_id.isoform_id
            BatchID = self.find_Next_BatchID(IsoformID, self.batch_id)
            if BatchID:
                self.batch_id = BatchID
                self.isoformbatch_id = self.str_IsoBatchID(IsoformID, BatchID)
                super(Isoform_Batch, self).save(*args, **kwargs)
        else:
            # confirms Batch_ID from IsoBatchID
            self.batch_id = str(self.isoformbatch_id).split(PEPBATCH_SEP)[-1]
            super(Isoform_Batch, self).save(*args, **kwargs)
        
# ================================================================================================
            
#=================================================================================================
class Modification_Dictionary(AuditModel):
    """
    Biologic Modification Dictionary
    """
#=================================================================================================

    ID_SEQUENCE = 'BioMod_Dictionary'
    ID_PREFIX = 'BMOD'
    ID_PAD = 6

    LIST_VIEW_FIELDS = {
        'biomod_id': 'BioMod ID',
        'biomod_code': 'BioMod Code',
        'mw': 'MW',
        'nseq': 'nSeq',
        'fasta': 'Sequence',
    }

    biomod_id = models.CharField(primary_key=True, max_length=15, verbose_name = "BioMod ID")
    biomod_code = models.CharField(max_length=25, verbose_name = "BioMod Code")

    biln = models.TextField(blank=True, verbose_name = "BILN") 
    helm = models.TextField(blank=True, verbose_name = "HELM") 
    fasta = models.TextField(blank=True, verbose_name = "FASTA") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    nseq = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nSeq")

    class Meta:
        app_label = 'dbiologic'
        db_table = 'biomod_dictionary'
        ordering=['biomod_code']
        indexes = [
            models.Index(name="biomod_id_idx",fields=['biomod_id']),
            models.Index(name="biomod_code_idx",fields=['biomod_code']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,BIOMODID, BIOMODCode=None, verbose=0):
    # Returns an instance by compound_id
        try:
            if BIOMODID:
                retInstance = cls.objects.get(biomod_id=BIOMODID)
            elif BIOMODCode:
                retInstance = cls.objects.get(biomod_code=BIOMODCode)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if BIOMODID:
                    logger.warning(f"[Biologic Modification Not Found] {BIOMODID} ")
                elif BIOMODCode :
                    logger.warning(f"[Biologic Modification Not Found] {BIOMODCode}  ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,BIOMODID, BIOMODCode=None, verbose=0):
    # Returns if an instance exists by compound_id
        if BIOMODID:
            retValue = cls.objects.filter(biomod_id=BIOMODID).exists()
        elif BIOMODCode:
            retValue = cls.objects.filter(biomod_code=BIOMODCode).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        _sync_sequence_fields(self)
        _compute_mw_nseq(self)
        if not self.biomod_id:
            self.biomod_id = self.next_id()
            if self.biomod_id:
                #self.save_batch() 
                super(Modification_Dictionary, self).save(*args, **kwargs)
        else:
            #self.save_batch()
            super(Modification_Dictionary, self).save(*args, **kwargs) 


#=================================================================================================
class Monomer(AuditModel):
    """
    Monomer library entry for BILN/HELM biologics.

    Mirrors the pyPept monomer SDF (symbol, m_type, m_Rgroups, ...) so that
    the monomer library can be stored in the database, edited by users, and
    read back by pyPept via ``spec_from_database``.
    """
#=================================================================================================

    ID_SEQUENCE = 'Monomer'
    ID_PREFIX = 'MON'
    ID_PAD = 5

    LIST_VIEW_FIELDS = {
        'monomer_id': 'Monomer ID',
        'symbol': 'Symbol',
        'monomer_name': 'Name',
        'monomer_type': 'Type',
        'monomer_subtype': 'SubType',
        'nat_analog': 'Natural Analog',
        'pdb_name': 'PDB',
    }

    DETAIL_VIEW_GROUPS = [
        ['monomer_id', 'symbol', 'monomer_abbr', 'monomer_name'],
        ['monomer_type', 'monomer_subtype', 'nat_analog', 'pdb_name'],
    ]

    monomer_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Monomer ID")
    symbol = models.CharField(max_length=30, unique=True, verbose_name = "Symbol")
    monomer_name = models.CharField(max_length=150, blank=True, verbose_name = "Name")
    monomer_abbr = models.CharField(max_length=30, blank=True, verbose_name = "Abbreviation")
    monomer_type = models.CharField(max_length=15, blank=True, verbose_name = "Type")
    monomer_subtype = models.CharField(max_length=20, blank=True, verbose_name = "SubType")
    nat_analog = models.CharField(max_length=10, blank=True, verbose_name = "Natural Analog")
    pdb_name = models.CharField(max_length=10, blank=True, verbose_name = "PDB Name")

    # Chemical structure: mol for querying/rendering, molfile keeps the exact
    # SDF mol block (including R-group annotations) for pyPept round-tripping.
    mol = rdm.MolField(blank=True, null=True, verbose_name = "Structure")
    molfile = models.TextField(blank=True, verbose_name = "Mol Block")
    
    rgroups = ArrayField(models.CharField(max_length=10, blank=True, null=True),
                         size=4, null=True, blank=True, verbose_name = "R-Groups")
    rgroup_idx = ArrayField(models.IntegerField(null=True, blank=True),
                            size=4, null=True, blank=True, verbose_name = "R-Group Idx")
    attachment_idx = ArrayField(models.IntegerField(null=True, blank=True),
                                size=4, null=True, blank=True, verbose_name = "Attachment Idx")

    class Meta:
        app_label = 'dbiologic'
        db_table = 'monomer'
        ordering = ['symbol']
        indexes = [
            models.Index(name="monomer_symbol_idx", fields=['symbol']),
            models.Index(name="monomer_type_idx", fields=['monomer_type']),
            GistIndex(name="monomer_mol_idx", fields=['mol']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.symbol} ({self.monomer_id})"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.symbol}"

    #------------------------------------------------
    @classmethod
    def get(cls, Symbol, verbose=0):
        try:
            retInstance = cls.objects.get(symbol=Symbol)
        except cls.DoesNotExist:
            if verbose:
                logger.warning(f"[Monomer Not Found] {Symbol} ")
            retInstance = None
        return retInstance

    #------------------------------------------------
    @classmethod
    def exists(cls, Symbol, verbose=0):
        return cls.objects.filter(symbol=Symbol).exists()

    #------------------------------------------------
    def set_molecule(self, molblock):
        """Parse a mol block and set both mol and molfile."""
        from rdkit import Chem
        mol = Chem.MolFromMolBlock(molblock) if molblock else None
        if mol is not None:
            self.mol = mol
            self.molfile = molblock
        return mol

    #------------------------------------------------
    def save(self, *args, **kwargs):
        # Keep the MolField in sync with the editable molfile.
        if self.molfile:
            from rdkit import Chem
            mol = Chem.MolFromMolBlock(self.molfile)
            if mol is not None:
                self.mol = mol
        if not self.monomer_id:
            self.monomer_id = self.next_id()
            if self.monomer_id:
                super(Monomer, self).save(*args, **kwargs)
        else:
            super(Monomer, self).save(*args, **kwargs)

    #------------------------------------------------
    @classmethod
    def spec_from_database(cls):
        """Return a pyPept-compatible monomer DataFrame read from the database.

        Equivalent to ``pyPept.sequence.get_monomer_info(sdf_path)``: indexed by
        ``symbol`` with a ``m_romol`` column holding RDKit molecules, plus the
        monomer property columns (m_Rgroups/m_RgroupIdx/m_attachmentPointIdx
        parsed into lists).
        """
        import pandas as pd
        from rdkit import Chem

        records = []
        for obj in cls.objects.filter(astatus__gte=0).order_by('symbol'):
            mol = None
            if obj.molfile:
                mol = Chem.MolFromMolBlock(obj.molfile)
            if mol is None and obj.mol is not None:
                mol = obj.mol

            records.append({
                'symbol': obj.symbol,
                'm_name': obj.monomer_name,
                'm_abbr': obj.monomer_abbr,
                'm_type': obj.monomer_type,
                'm_subtype': obj.monomer_subtype,
                'm_Rgroups': list(obj.rgroups or []) or [None, None, None, None],
                'm_RgroupIdx': list(obj.rgroup_idx or []) or [None, None, None, None],
                'm_attachmentPointIdx': list(obj.attachment_idx or []) or [None, None, None, None],
                'natAnalog': obj.nat_analog,
                'pdbName': obj.pdb_name,
                'm_romol': mol,
            })

        df = pd.DataFrame(records)
        if df.empty:
            return df
        return df.set_index('symbol')


def spec_from_database():
    """Module-level alias for Monomer.spec_from_database()."""
    return Monomer.spec_from_database()


#=================================================================================================
class Biologic_Modification(models.Model):
    """
    Ordered modification link between a Biologic_Isoform and a
    Modification_Dictionary entry, split by terminal position.
    """
#=================================================================================================
    MOD_POSITION_CHOICES = [
        ('N', 'N-Terminal'),
        ('C', 'C-Terminal'),
        ('I', 'Internal'),
    ]

    isoform_id = models.ForeignKey(Biologic_Isoform, null=False, blank=False, verbose_name = "Isoform ID", on_delete=models.CASCADE,
        db_column="isoform_id", related_name="modification_links")
    biomod_id = models.ForeignKey(Modification_Dictionary, null=False, blank=False, verbose_name = "Modification", on_delete=models.CASCADE,
        db_column="biomod_id", related_name="isoform_links")
    mod_position = models.CharField(max_length=1, choices=MOD_POSITION_CHOICES, default='I', verbose_name = "Position")
    mod_order = models.PositiveSmallIntegerField(default=0, verbose_name = "Order")

    class Meta:
        app_label = 'dbiologic'
        db_table = 'biologic_modification'
        ordering = ['mod_position', 'mod_order']
        unique_together = ('isoform_id', 'mod_position', 'mod_order')

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.isoform_id} <-[{self.mod_position}:{self.mod_order}]- {self.biomod_id}"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.biomod_id} ({self.get_mod_position_display()} {self.mod_order})"

    #------------------------------------------------
    def save(self, *args, **kwargs):
        super(Biologic_Modification, self).save(*args, **kwargs)
        # Rebuild the parent isoform's sequence whenever a modification link
        # is added or changed. Best-effort: a rebuild failure must not roll
        # back the saved modification link.
        try:
            self.isoform_id.rebuild_sequence_from_modifications()
        except Exception:
            logger.exception('Failed to rebuild isoform sequence after modification save')

    #------------------------------------------------
    def delete(self, *args, **kwargs):
        isoform = self.isoform_id
        super(Biologic_Modification, self).delete(*args, **kwargs)
        try:
            isoform.rebuild_sequence_from_modifications()
        except Exception:
            logger.exception('Failed to rebuild isoform sequence after modification removal')

