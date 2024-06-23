import re
from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError

#from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Chemical Structures and Samples 
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Library(AuditModel):
    """
    List of Sample Library
    """
#=================================================================================================
    Choice_Dictionary = {
        'library_class':'Library_Class',
    }
    ID_SEQUENCE = 'ChemLibrary'
    ID_PREFIX = 'LIB'
    ID_PAD = 5

    library_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Library ID")
    library_name = models.CharField(max_length=50, unique=True, verbose_name = "Library Name")
    chem_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="chem_class", related_name="%(class)s_chemclass")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    class Meta:
        app_label = 'dchem'
        db_table = 'library'
        ordering=['library_name']
        indexes = [
            models.Index(name="lib_lname_idx", fields=['library_name']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.library_name} ({self.library_id}) {self.source}"

    #------------------------------------------------
    @classmethod
    def get(cls,LibraryID=None,LibraryName=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if LibraryID:
                retInstance = cls.objects.get(library_id=LibraryID)
            elif LibraryName:
                retInstance = cls.objects.get(library_name=LibraryName)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if LibraryID:
                    print(f"[Library Not Found] {LibraryID} ")
                elif LibraryName:
                    print(f"[Library Not Found] {LibraryName} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,LibraryID=None,LibraryName=None,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        if LibraryID:
            retValue = cls.objects.filter(library_id=LibraryID).exists()
        elif LibraryName:
            retValue = cls.objects.filter(library_name=LibraryName).exists()
        else:
            retValue = False
        return(retValue)


    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.library_id:
            self.library_id = self.next_id()
            if self.library_id: 
                super(Library, self).save(*args, **kwargs)
        else:
            super(Library, self).save(*args, **kwargs) 


#=================================================================================================
class Chem_Structure(AuditModel):
    """
    List of ChemStructure 
    """
#=================================================================================================
    Choice_Dictionary = {
        'chem_class':'Chem_Class',
    }

    ID_SEQUENCE = 'ChemStructure'
    ID_PREFIX = 'CS'
    ID_PAD = 7

    structure_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Structure ID")
    structure_name = models.CharField(max_length=50, unique=True, verbose_name = "Structure Name")
    chem_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
        db_column="chem_class", related_name="%(class)s_chemclass")

    smol = models.MolField(blank=True, null=True, verbose_name = "MOL")	
    torsionbv = models.BfpField(null=True)	
    ffp2 = models.BfpField(null=True, verbose_name = "FFP2")
    mfp2 = models.BfpField(null=True, verbose_name = "MFP2")

    mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
    mw = models.FloatField(default=0, blank=True, verbose_name ="MW")

    class Meta:
        app_label = 'dchem'
        db_table = 'chemstructure'
        ordering=['structure_name']
        indexes = [
            models.Index(name="cstruct_dname_idx", fields=['structure_name']),
            GistIndex(name="cstruct_smol_idx",fields=['smol']),
            GistIndex(name="cstruct_ffp2_idx",fields=['ffp2']),
            GistIndex(name="cstruct_mfp2_idx",fields=['mfp2'])
        ]

    #------------------------------------------------
    # def __str__(self) -> str:
    #     return f"{self.drug_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.structure_name} ({self.structure_id})"

    #------------------------------------------------
    @classmethod
    def get(cls,StructureID=None,StructureName=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if StructureID:
                retInstance = cls.objects.get(structure_id=StructureID)
            elif StructureName:
                retInstance = cls.objects.get(structure_name=StructureName)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if StructureID:
                    print(f"[Structure Not Found] {StructureID} ")
                elif StructureName:
                    print(f"[Structure Not Found] {StructureName} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,StructureID=None,StructureName=None,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        if StructureID:
            retValue = cls.objects.filter(structure_id=StructureID).exists()
        elif StructureName:
            retValue = cls.objects.filter(structure_name=StructureName).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    @classmethod
    def update_all_fp(cls):
        # Require RDKit-PostgreSQL
        cls.objects.update(ffp2=FEATMORGANBV_FP('smol'),
                           mfp2=MORGANBV_FP('smol'), 
                           torsionbv=TORSIONBV_FP('smol')
                           )

    #------------------------------------------------
    @classmethod
    def smiles2mol(cls,Smiles,verbose=0):
        try:
            xmol = Chem.MolFromSiles(Smiles)
        except:
            xmol = None
            if verbose:
                print(f"[Invalid SMILES] {Smiles} ")
        return(xmol)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.structure_id:
            self.structure_id = self.next_id()
            if self.structure_id: 
                super(Chem_Structure, self).save(*args, **kwargs)
                # self.__dict__.update(ffp2=FEATMORGANBV_FP('smol'), mfp2=MORGANBV_FP('smol'), torsionbv=TORSIONBV_FP('smol'))

                # Require RDKit-PostgreSQL
                # ChemStructure.objects.filter(structure_id=self.structure_id).update(
                #     ffp2=FEATMORGANBV_FP('smol'), 
                #     mfp2=MORGANBV_FP('smol'), 
                #     torsionbv=TORSIONBV_FP('smol')
                #     )
        else:
            # Require RDKit-PostgreSQL
            # self.__dict__.update(
            #     ffp2=FEATMORGANBV_FP('smol'), 
            #     mfp2=MORGANBV_FP('smol'), 
            #     torsionbv=TORSIONBV_FP('smol')
            #     )
            super(Chem_Structure, self).save(*args, **kwargs) 

#=================================================================================================
class Sample(AuditModel):
    """
    List of ChemStructure 
    """
#=================================================================================================
    ID_SEQUENCE = 'Sample'
    ID_PREFIX = 'S'
    ID_PAD = 7

    sample_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Sample ID")
    sample_name = models.CharField(max_length=50, unique=True, verbose_name = "Sample Name")
    library_id = models.ForeignKey(Library, null=False, blank=False, verbose_name = "Library ID", on_delete=models.DO_NOTHING,
        db_column="library_id", related_name="%(class)s_library_id")

    sample_code = models.CharField(max_length=50, unique=True, verbose_name = "Sample Code")
    structure_id = models.ForeignKey(Chem_Structure, null=False, blank=False, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
        db_column="structure_id", related_name="%(class)s_structure_id")
    other_ids = models.CharField(max_length=50, unique=True, verbose_name = "Other IDs")

    chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
        db_column="chemist", related_name="%(class)s_chemist")

    class Meta:
        app_label = 'dchem'
        db_table = 'sample'
        ordering=['library_id','sample_name']
        indexes = [
            models.Index(name="smp_sname_idx", fields=['sample_name']),
            models.Index(name="smp_scode_idx", fields=['sample_code']),
        ]

#=================================================================================================
class Sample_Batch(AuditModel):
    """
    List of Sample Batches 
    """
#=================================================================================================

    samplebatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "SampleBatch ID")
    sample_id = models.ForeignKey(Sample, null=False, blank=False, verbose_name = "Sample ID", on_delete=models.DO_NOTHING,
        db_column="sample_id", related_name="%(class)s_sample_id")
    previous_batch_id= models.CharField(max_length=20, blank=True, verbose_name = "Previous SampleBatch ID")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[AlphaNumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")

    salt = models.CharField(max_length=500, blank=True, verbose_name = "Salt")
    mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
    mw = models.FloatField(default=0, blank=True, verbose_name ="MW")

    quality_source = models.CharField(max_length=150, blank=True, verbose_name = "QC Source")
    qc_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "QC status", on_delete=models.DO_NOTHING,
        db_column="qc_status", related_name="%(class)s_qc")
    qc_record = models.CharField(max_length=150, blank=True, verbose_name = "QC Records")
    stock_date = models.DateField(null=True, blank=True, verbose_name = "Stock Date") 
    stock_level = models.CharField(max_length=20, blank=True, verbose_name = "Stock Levels") 
    chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
        db_column="chemist", related_name="%(class)s_chemist")

    #------------------------------------------------
    class Meta:
        app_label = 'dchem'
        db_table = 'samplebatch'
        ordering=['samplebatch_id']
        indexes = [
            models.Index(name="smpbatch_samplebatch_idx",fields=['sample_id','batch_id']),
            models.Index(name="smpbatch_qc_idx",fields=['qc_status']),
            models.Index(name="smpbatch_sdate_idx",fields=['stock_date']),
            models.Index(name="smpbatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.samplebatch_id}"

    #------------------------------------------------
    @classmethod
    # Formats BatchNo:int -> BatchID:str 
    def str_BatchID(self,BatchNo:int) -> str:
        return(f"{BatchNo:02d}")
    #------------------------------------------------
    @classmethod
    def str_SampleBatchID(self,SampleID:str,BatchID:str) -> str:
        return(f"{SampleID}{SAMPLEBATCH_SEP}{BatchID}")

    #------------------------------------------------
    def find_Next_BatchID(self, SampleID:str, BatchID:str=None) -> str:
        # Check for given BatchID    
        if BatchID:
            # Clean up BatchID - remove non alphanumeric character and make uppercase
            BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

            # Clean up BatchID - reformat numbers
            if BatchID.isnumeric():
                BatchID = self.str_BatchID(int(BatchID))

            next_SampleBatch = self.str_SampleBatchID(SampleID,BatchID)
            if not self.exists(next_SampleBatch):
                return(BatchID)

        # Find new BatchID    
        next_BatchNo = 1
        next_SampleBatch = self.str_SampleBatchID(SampleID,self.str_BatchID(next_BatchNo))
        while self.exists(next_SampleBatch):
            next_BatchNo = next_BatchNo + 1
            next_SampleBatch = self.str_SampleBatchID(SampleID,self.str_BatchID(next_BatchNo))
        return(self.str_BatchID(next_BatchNo))    

    #------------------------------------------------
    @classmethod
    def get(cls,SampleBatchID,verbose=0):
    # Returns an instance if found by samplebatch_id
        try:
            retInstance = cls.objects.get(samplebatch_id=SampleBatchID)
        except:
            if verbose:
                print(f"[SampleBatch Not Found] {SampleBatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,SampleBatchID,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(samplebatch_id=SampleBatchID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.samplebatch_id: 
            # creates new SampleBatchID
            SampleID = self.sample_id.sample_id
            BatchID = self.find_Next_BatchID(SampleID,self.batch_id)
            if BatchID:
                self.batch_id = BatchID
                self.samplebatch_id = self.str_SampleBatchID(SampleID,BatchID)
                super(Sample_Batch,self).save(*args, **kwargs)
        else:
            # confirms Batch_ID from SampleBatchID
            self.batch_id = str(self.samplebatch_id).replace(str(self.sample_id.sample_id),"").split(SAMPLEBATCH_SEP)[1]
            super(Sample_Batch,self).save(*args, **kwargs)
            #print(f"[SampleBatch.save]: {self.samplebatch_id}")
