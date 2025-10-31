import re

#from django_rdkit import models
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from dorganism.models import Taxonomy 
from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Peptide Application Model

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
        'seq':'Peptide Sequence',
        # 'peptide_name':'Peptide Name',
        # 'peptide_type':'Peptide Type',
        # 'peptide_panel':'Panel',
        # 'peptide_notes':'Notes',
        # 'peptide_code':'Peptide Code',
        # 'peptide_identification':'Identification',
        # 'peptide_origin':'Origin',
        # 'source':"Source",
        # 'source_code':"Source Code",
        # 'reference': "Reference",
        #'tax_id':{'Tax-ID': {'tax_id':URL_LINKS['tax_id']}}, 
    }

    CARDS_FIELDS= {
        "source" : "Source",
        "peptide_type": "Type",
    #     "res_property": "Phenotype",
    #     "gen_property": "Genotype",
    }

    VIEW_GROUPS = [
        ['peptide_name', 'peptide_type', 'peptide_panel', 'peptide_notes'],
        ['seq','bilm','helm'],
        ['peptide_origin','source', 'source_code','reference'],
        ['mta_status','mta_notes','mta_document','biologist'],
    ]

    DICTIONARY_FIELDS = {
        'mta_status':'License_Status',
        'peptide_type':'Peptide_Type',
        'peptide_panel':'Peptide_Panel',
    }

    peptide_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Peptide ID") 
    peptide_name= models.CharField(max_length=200, blank=True, verbose_name = "Peptide Name") 
    peptide_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Peptide Notes")
    #peptide_code= models.CharField(max_length=30, blank=True, verbose_name = "Peptide Code")
    peptide_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    peptide_type=models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Peptide Type", on_delete=models.DO_NOTHING,
        db_column="peptide_type", related_name="%(class)s_peptide_type")
    #peptide_identification = models.CharField(max_length=512, blank=True, verbose_name = "Peptide Identification")
    peptide_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Peptide")
    bilm= models.CharField(max_length=200, blank=True, verbose_name = "BILM") 
    helm= models.CharField(max_length=200, blank=True, verbose_name = "HELM") 
    seq= models.CharField(max_length=200, blank=True, verbose_name = "Seq") 

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
            models.Index(name="pep_type_idx", fields=['peptide_type']),
            models.Index(name="pep_panel_idx", fields=['peptide_panel']),
            #models.Index(name="pep_source_idx", fields=['source']),
            # models.Index(name="org_taxid_idx", fields=['tax_id']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.peptide_id} ({self.peptide_name})"

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
                retInstance = cls.objects.get(peptide_ID=PeptideID)
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
        return(f"{PeptidClass}{PEPTIDE_SEP}{PeptideNo:04d}")


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
            self.peptide_id = self.find_Next_PeptideID()
            if self.peptide_id: 
                super(Peptide, self).save(*args, **kwargs)
        else:
            super(Peptide, self).save(*args, **kwargs) 

#=================================================================================================
class Peptide_Batch(AuditModel):
    """
    Peptide Batch Collection
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
       'Group1': ["batch_id", "batch_notes", "previous_batch_id", "passage_number", "qc_status", "batch_quality", "quality_source", "stock_date", "stock_level", "biologist" ]
       }



    pepbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "PepBatch ID")
    peptide_id = models.ForeignKey(Peptide, null=False, blank=False, verbose_name = "Peptide ID", on_delete=models.DO_NOTHING,
        db_column="peptide_id", related_name="%(class)s_peptide_id")
    #previous_batch_id= models.CharField(max_length=20, blank=True, verbose_name = "Previous peptideBatch ID")
    #passage_number= models.CharField(max_length=20, blank=True, verbose_name = "Passage Number")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    source_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Source Type", on_delete=models.DO_NOTHING,
        db_column="source_type", related_name="%(class)s_source_type")

    peptide_tags = 	models.CharField(max_length=80, blank=True, verbose_name = "Tags")
    # Expression system	
    # Expression cell ID	
    # Vector ID	Free form ID 1	Free form ID 2	
  
    # batch_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Quality", on_delete=models.DO_NOTHING,
    #     db_column="batch_quality", related_name="%(class)s_batchquality")

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
        db_table = 'pepbatch'
        ordering=['pepbatch_id']
        indexes = [
            models.Index(name="pepbatch_pepbatch_idx",fields=['peptide_id','batch_id']),
            models.Index(name="pepbatch_qc_idx",fields=['qc_status']),
            models.Index(name="pepbatch_sdate_idx",fields=['stock_date']),
            models.Index(name="pepbatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.pepbatch_id}"

    #------------------------------------------------
    @classmethod
    # Formats BatchNo:int -> BatchID:str 
    def str_BatchID(self,BatchNo:int) -> str:
        return(f"{BatchNo:02d}")
    #------------------------------------------------
    @classmethod
    # Formats PeptideID:str,BatchID:str -> PepBatchID:str
    def str_PepBatchID(self,PeptideID:str,BatchID:str) -> str:
        return(f"{PeptideID}{PEPTIDE_SEP}{BatchID}")

    #------------------------------------------------
    def find_Next_BatchID(self, PeptideID:str, BatchID:str=None) -> str:
        # Check for given BatchID    
        if BatchID:
            # Clean up BatchID - remove non alphanumeric character and make uppercase
            BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

            # Clean up BatchID - reformat numbers
            if BatchID.isnumeric():
                BatchID = self.str_BatchID(int(BatchID))

            next_PepBatch = self.str_PepBatchID(PeptideID,BatchID)
            if not self.exists(next_PepBatch):
                return(BatchID)

        # Find new BatchID    
        next_BatchNo = 1
        next_PepBatch = self.str_PepBatchID(PeptideID,self.str_BatchID(next_BatchNo))
        while self.exists(next_PepBatch):
            next_BatchNo = next_BatchNo + 1
            next_PepBatch = self.str_PepBatchID(PeptideID,self.str_BatchID(next_BatchNo))
        return(self.str_BatchID(next_BatchNo))    

    #------------------------------------------------
    @classmethod
    def get(cls,PepBatchID,verbose=0):
    # Returns an instance if found by pepbatch_id
        try:
            retInstance = cls.objects.get(pepbatch_id=PepBatchID)
        except:
            if verbose:
                print(f"[PepBatch Not Found] {PepBatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,PepBatchID,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(pepbatch_id=PepBatchID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.pepbatch_id: 
            # creates new PepBatchID
            PeptideID = self.peptide_id.peptide_id
            BatchID = self.find_Next_BatchID(PeptideID,self.batch_id)
            if BatchID:
                self.batch_id = BatchID
                self.pepbatch_id = self.str_PepBatchID(PeptideID,BatchID)
                super(Peptide_Batch,self).save(*args, **kwargs)
        else:
            # confirms Batch_ID from PepBatchID
            self.batch_id = str(self.pepbatch_id).replace(str(self.peptide_id.peptide_id),"").split(ORGBATCH_SEP)[1]
            super(Peptide_Batch,self).save(*args, **kwargs)
            #print(f"[PepBatch.save]: {self.pepbatch_id}")
        
# ================================================================================================


# Create your models here.
