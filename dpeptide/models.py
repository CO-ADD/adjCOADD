import re

#from django_rdkit import models
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from dorganism.models import Taxonomy 
from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------------------------
# Peptide Application Model
# Peptide
#  -> Peptide_Isoform  <- Modification_Dictionary (via Isoform_Modification)
#     -> Isoform_Batch
#
# AminoAcid (dchem.models) - List of Amino Acids for Peptide Sequences, including modified and non-natural amino acids
# Modification_Dictionary - List of Peptide Modifications/Tags, including modified and non-natural amino acids
#

#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Peptide(AuditModel):
    """
    Main class of Peptides
    
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
#       'organism_name':{"VerboseName":'Peptide Name','Updatable':False}
        'peptide_id':{'Peptide ID': {'peptide_id':URL_LINKS['peptide_id']}}, 
        'peptide_name':'Name',
        'peptide_class':'Type',
        # 'peptide_panel':'Panel',
        'peptide_notes':'Notes',
        # 'peptide_code':'Peptide Code',
        # 'peptide_identification':'Identification',
        # 'peptide_origin':'Origin',
        'seq':'Peptide Sequence',
        # 'source':"Source",
        # 'source_code':"Source Code",
        # 'reference': "Reference",
        #'tax_id':{'Tax-ID': {'tax_id':URL_LINKS['tax_id']}}, 
    }

    CARDS_FIELDS= {
        "source" : "Source",
        "peptide_class": "Class",
    #     "res_property": "Phenotype",
    #     "gen_property": "Genotype",
    }

    DETAIL_VIEW_GROUPS = [
        ['peptide_name', 'peptide_class', 'peptide_panel', 'peptide_notes'],
        ['seq','bilm','helm'],
        ['peptide_origin','source', 'source_code','reference'],
        ['mta_status','mta_notes','mta_document','biologist'],
    ]

    DICTIONARY_FIELDS = {
        'mta_status':'License_Status',
        'peptide_class':'Peptide_Class',
        'peptide_panel':'Peptide_Panel',
    }

    peptide_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Peptide ID") 
    peptide_name= models.CharField(max_length=200, blank=True, verbose_name = "Peptide Name") 
    peptide_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Peptide Notes")
    #peptide_code= models.CharField(max_length=30, blank=True, verbose_name = "Peptide Code")
    peptide_panel=ArrayField(models.CharField(max_length=100, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    peptide_class=models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="peptide_class", related_name="%(class)s_peptide_class")
    #peptide_identification = models.CharField(max_length=512, blank=True, verbose_name = "Peptide Identification")
    peptide_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Peptide")

    bilm= models.CharField(max_length=200, blank=True, verbose_name = "BILM") 
    helm= models.CharField(max_length=200, blank=True, verbose_name = "HELM") 
    seq= models.CharField(max_length=200, blank=True, verbose_name = "Seq") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    nseq = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nSeq")

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
        db_table = "peptide_doc", related_name="%(class)s_document")
    
    #tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")


#------------------------------------------------

    class Meta:
        app_label = 'dpeptide'
        db_table = 'peptide'
        #ordering=['peptide_id']
        indexes = [
            models.Index(name="pep_bilm_idx", fields=['bilm']),
            #models.Index(name="pep_code_idx", fields=['peptide_code']),
            models.Index(name="pep_type_idx", fields=['peptide_class']),
            GinIndex(name="pep_panel_idx", fields=['peptide_panel']),
            #models.Index(name="pep_source_idx", fields=['source']),
            # models.Index(name="org_taxid_idx", fields=['tax_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.peptide_id} ({self.peptide_name})"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.peptide_id} {self.peptide_name}".strip()

    #------------------------------------------------
    @classmethod
    def exists(cls,PeptideID=None,PeptideName=None,verbose=0):
        if PeptideID:
            # Returns if an instance exists by peptide_id
            return cls.objects.filter(peptide_id=PeptideID).exists()
        elif PeptideName:
            # Returns if an instance exists by peptide_line
            return cls.objects.filter(peptide_name=PeptideName).exists()

    #------------------------------------------------
    @classmethod
    def get(cls,PeptideID=None,PeptideName=None,verbose=0):
        if PeptideID:
            # Returns an instance by peptide_id
            try:
                retInstance = cls.objects.get(peptide_id=PeptideID)
            except:
                if verbose:
                    print(f"[PeptideID Not Found] {PeptideID} ")
                retInstance = None
        elif PeptideName:
            # Returns an instance by peptide_id
            try:
                retInstance = cls.objects.get(peptide_name=PeptideName)
            except:
                if verbose:
                    print(f"[Peptide Name Not Found] {PeptideName} ")
                retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def str_PeptideID(cls,PeptidClass,PeptideNo) -> str:
    #
    # Output:   Peptide_ID as string like MAB_0001, PEP_0001, NB_0001 
    #
        return(f"{PeptidClass}{PEPID_SEP}{PeptideNo:04d}")


    #------------------------------------------------
    @classmethod
    def find_Next_PeptideID(cls,PeptidClass, PeptideClassType = PEPTIDE_CLASSES) -> str:
        if PeptidClass in PeptideClassType:
            Peptide_IDSq=Sequence(PeptidClass)
            Peptide_nextID = next(Peptide_IDSq)
            Peptide_strID = cls.str_PeptideID(PeptidClass,Peptide_nextID)
            while cls.exists(Peptide_strID):
                Peptide_nextID = next(Peptide_IDSq)
                Peptide_strID = cls.str_PeptideID(PeptidClass,Peptide_nextID)
            return(Peptide_strID)    
        else:
            return(None)

     #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.peptide_id:
            _class = str(self.peptide_class.dict_value) if self.peptide_class else 'PEP'
            self.peptide_id = self.find_Next_PeptideID(_class)
            if self.peptide_id:
                super(Peptide, self).save(*args, **kwargs)
        else:
            super(Peptide, self).save(*args, **kwargs) 


#=================================================================================================
class Peptide_Isoform(AuditModel):
    """
    Isoform - a modified variant of a Peptide
    
    """
#=================================================================================================
    LIST_VIEW_FIELDS = {
        'isoform_id': {'Isoform ID': {'isoform_id': URL_LINKS['isoform_id']}},
        'isoform_name': 'Isoform Name',
        'peptide_id': 'Peptide',
        #'seq': 'Sequence',
        'mw': 'MW',
        'biologist': 'Biologist',
    }

    DETAIL_VIEW_GROUPS = [
        ['isoform_id', 'isoform_name', 'isoform_notes', 'biologist'],
        ['peptide_id', 'mw', 'nseq'],
        ['seq', 'bilm', 'helm', ],
    ]

    DICTIONARY_FIELDS = {}

    isoform_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Isoform ID") 
    isoform_name= models.CharField(max_length=200, blank=True, verbose_name = "Isoform Name") 
    isoform_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Isoform Notes")
    
    peptide_id = models.ForeignKey(Peptide, null=False, blank=False, verbose_name = "Peptide ID", on_delete=models.DO_NOTHING,
        db_column="peptide_id", related_name="isoforms")

    peptide_modifications = models.ManyToManyField(
        'Modification_Dictionary',
        through='Isoform_Modification',
        through_fields=('isoform_id', 'pepmod_id'),
        blank=True,
        related_name='isoforms',
        verbose_name = "Modifications",
    )

    bilm= models.CharField(max_length=200, blank=True, verbose_name = "BILM") 
    helm= models.CharField(max_length=200, blank=True, verbose_name = "HELM") 
    seq= models.CharField(max_length=200, blank=True, verbose_name = "Seq") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    nseq = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nSeq")

    assoc_documents = models.ManyToManyField(Document,verbose_name = "Documents", blank=True,
        db_table = "isoform_doc", related_name="%(class)s_document")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    class Meta:
        app_label = 'dpeptide'
        db_table = 'peptide_isoform'
        ordering=['peptide_id','isoform_id']
        indexes = [
            models.Index(name="pepiso_bilm_idx", fields=['bilm']),
            models.Index(name="pepiso_pep_idx", fields=['peptide_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.isoform_id} ({self.isoform_name})"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.isoform_id} {self.isoform_name}".strip()

    #------------------------------------------------
    @classmethod
    def str_IsoformID(cls, PeptideID, IsoformNo) -> str:
        return(f"{PeptideID}{PEPISO_SEP}I{IsoformNo:02d}")

    #------------------------------------------------
    def find_Next_IsoformID(self) -> str:
        peptide_id = self.peptide_id.peptide_id
        next_no = 1
        next_iso = self.str_IsoformID(peptide_id, next_no)
        while self.exists(next_iso):
            next_no += 1
            next_iso = self.str_IsoformID(peptide_id, next_no)
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
        if not self.isoform_id:
            self.isoform_id = self.find_Next_IsoformID()
            if self.isoform_id:
                super(Peptide_Isoform, self).save(*args, **kwargs)
        else:
            super(Peptide_Isoform, self).save(*args, **kwargs)

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
    Peptide Isoform Batch Collection
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
        "stock_date":"Stock Date",
        "stock_level":"Stock Levels",
        "biologist":"Biologist"
    }

    DICTIONARY_FIELDS = {
        'qc_status':'QC_Status',
        'source_type' : 'Peptide_Source'
    }
    
    FORM_GROUPS = {
       'Group1': ["batch_id", "batch_notes", "qc_status", "quality_source", "stock_date", "stock_level", "biologist", "source", "source_code", "source_type" ]
       }

    isoformbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "IsoFormBatch ID")
    isoform_id = models.ForeignKey(Peptide_Isoform, null=False, blank=False, verbose_name = "Isoform ID", on_delete=models.DO_NOTHING,
        db_column="isoform_id", related_name="batches")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    source_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Source Type", on_delete=models.DO_NOTHING,
        db_column="source_type", related_name="%(class)s_source_type")

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
        app_label = 'dpeptide'
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
    Peptide Modification Dictionary
    """
#=================================================================================================

    ID_SEQUENCE = 'PepMod_Dictionary'
    ID_PREFIX = 'PMOD'
    ID_PAD = 6

    LIST_VIEW_FIELDS = {
        'pepmod_id': 'PepMod ID',
        'pepmod_code': 'PepMod Code',
        'smiles': 'SMILES',
        'mf': 'Formula',
        'mw': 'MW',
    }

    pepmod_id = models.CharField(primary_key=True, max_length=15, verbose_name = "PepMod ID")
    pepmod_code = models.CharField(max_length=25, verbose_name = "PepMod Code")

    bilm= models.CharField(max_length=200, blank=True, verbose_name = "BILM") 
    helm= models.CharField(max_length=200, blank=True, verbose_name = "HELM") 
    seq= models.CharField(max_length=200, blank=True, verbose_name = "Seq") 
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")

    class Meta:
        app_label = 'dpeptide'
        db_table = 'pepmod_dictionary'
        ordering=['pepmod_code']
        indexes = [
            models.Index(name="pepmod_id_idx",fields=['pepmod_id']),
            models.Index(name="pepmod_code_idx",fields=['pepmod_code']),
        ]

    #------------------------------------------------
    @classmethod
    def get(cls,PEPMODID, PEPMODCode=None, verbose=0):
    # Returns an instance by compound_id
        try:
            if PEPMODID:
                retInstance = cls.objects.get(pepmod_id=PEPMODID)
            elif PEPMODCode:
                retInstance = cls.objects.get(pepmod_code=PEPMODCode)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if PEPMODID:
                    logger.warning(f"[Peptide Modification Not Found] {PEPMODID} ")
                elif PEPMODCode :
                    logger.warning(f"[Peptide Modification Not Found] {PEPMODCode}  ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,PEPMODID, PEPMODCode=None, verbose=0):
    # Returns if an instance exists by compound_id
        if PEPMODID:
            retValue = cls.objects.filter(pepmod_id=PEPMODID).exists()
        elif PEPMODCode:
            retValue = cls.objects.filter(pepmod_code=PEPMODCode).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.pepmod_id:
            self.pepmod_id = self.next_id()
            if self.pepmod_id:
                #self.save_batch() 
                super(Modification_Dictionary, self).save(*args, **kwargs)
        else:
            #self.save_batch()
            super(Modification_Dictionary, self).save(*args, **kwargs) 


#=================================================================================================
class Isoform_Modification(models.Model):
    """
    Ordered modification link between a Peptide_Isoform and a
    Modification_Dictionary entry, split by terminal position.
    """
#=================================================================================================
    MOD_POSITION_CHOICES = [
        ('N', 'N-Terminal'),
        ('C', 'C-Terminal'),
        ('I', 'Internal'),
    ]

    isoform_id = models.ForeignKey(Peptide_Isoform, null=False, blank=False, verbose_name = "Isoform ID", on_delete=models.CASCADE,
        db_column="isoform_id", related_name="modification_links")
    pepmod_id = models.ForeignKey(Modification_Dictionary, null=False, blank=False, verbose_name = "Modification", on_delete=models.CASCADE,
        db_column="pepmod_id", related_name="isoform_links")
    mod_position = models.CharField(max_length=1, choices=MOD_POSITION_CHOICES, default='I', verbose_name = "Position")
    mod_order = models.PositiveSmallIntegerField(default=0, verbose_name = "Order")

    class Meta:
        app_label = 'dpeptide'
        db_table = 'isoform_modification'
        ordering = ['mod_position', 'mod_order']
        unique_together = ('isoform_id', 'mod_position', 'mod_order')

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.isoform_id} <-[{self.mod_position}:{self.mod_order}]- {self.pepmod_id}"

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.pepmod_id} ({self.get_mod_position_display()} {self.mod_order})"

