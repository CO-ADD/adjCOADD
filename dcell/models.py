from django.db import models
import re
from model_utils import Choices
from sequences import Sequence
#from django_rdkit import models
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from django.core.validators import RegexValidator
from dorganism.models import Taxonomy 

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Cell Application Model

#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Cell(AuditModel):
    """
    Main class of Cells in Isolate Collection
    
    """
#=================================================================================================
    HEADER_FIELDS = {
#       'organism_name':{"VerboseName":'Cell Name','Updatable':False}
        'cell_id':{'Cell ID': {'cell_id':LinkList['cell_id']}}, 
        'cell_names':'Cell Names',
        'cell_line':'Cell Line',
        'cell_type':'Cell Type',
        'cell_panel':'Panel',
        'cell_notes':'Notes',
        'cell_code':'Cell Code',
        'cell_idenfication':'Identification',
        'cell_origin':'Origin',
        'source':"Source",
        'source_code':"Source Code",
        'reference': "Reference",
        #'tax_id':{'Tax-ID': {'tax_id':LinkList['tax_id']}}, 
    }

    CARDS_FIELDS= {
        "source" : "Source",
        "cell_type": "Type",
    #     "res_property": "Phenotype",
    #     "gen_property": "Genotype",
    }

    FORM_GROUPS={
       'Group1': ["cell_line", "cell_type", "cell_names", "cell_panel", "cell_origin", "cell_notes"],
       'Group2': ['source', 'source_code','reference'],
       'Group3': ['mta_status','mta_notes','mta_document','biologist'],
       'Group4': ['collect_tissue', 'patient_diagnosis', 'patient']
    }

    Choice_Dictionary = {
        'mta_status':'License_Status',
        'cell_type':'Cell_Type',
        'cell_panel':'Cell_Panel',
    }

    cell_id = models.CharField(primary_key=True, max_length=15, verbose_name = "Cell ID") 
    cell_names= models.CharField(max_length=200, blank=True, verbose_name = "Cell Names") 
    cell_line= models.CharField(max_length=200, blank=True, verbose_name = "Cell Line") 
    cell_notes= models.CharField(max_length=1024, blank=True, verbose_name = "Cell Notes")
    cell_code= models.CharField(max_length=30, blank=True, verbose_name = "Cell Code")
    cell_panel=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)
    cell_type=ArrayField(models.CharField(max_length=100, null=True, blank=True), size=20, verbose_name = "Type", null=True, blank=True)
    cell_identification = models.CharField(max_length=512, blank=True, verbose_name = "Cell Identification")
    cell_origin = models.CharField(max_length=512, blank=True, verbose_name = "Origin of Cell")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
    reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

    organism_name= models.ForeignKey(Taxonomy, null=False, blank=False, verbose_name = "Organism Name", on_delete=models.DO_NOTHING, 
        db_column="organism_name", related_name="%(class)s_organism_name")

    mta_status = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MTA Status", on_delete=models.DO_NOTHING,
        db_column="mta_status", related_name="%(class)s_mta")
    mta_notes = models.CharField(max_length=150, blank=True, verbose_name = "MTA Notes")
    mta_document = models.CharField(max_length=150, blank=True, verbose_name = "MTA Document")
    mta_notes = models.CharField(max_length=512, blank=True, verbose_name = "MTA Notes")

    collect_tissue = models.CharField(max_length=120, blank=True, verbose_name = "From Tissue/Organ")
    patient_diagnosis = models.CharField(max_length=120, blank=True, verbose_name = "Patient Diagnosis")
    patient = models.CharField(max_length=20, blank=True, verbose_name = "Patient Info")

    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    assoc_documents = models.ManyToManyField(Document,verbose_name = "Documents", blank=True,
        db_table = "cell_doc", related_name="%(class)s_document")
    
    #tax_id = models.IntegerField(default=0, verbose_name = "NCBI Tax ID")


#------------------------------------------------


    class Meta:
        app_label = 'dcell'
        db_table = 'cell'
        #ordering=['cell_id']
        indexes = [
            models.Index(name="cell_line_idx", fields=['cell_line']),
            models.Index(name="cell_cellcode_idx", fields=['cell_code']),
            models.Index(name="cell_celltype_idx", fields=['cell_type']),
            models.Index(name="cell_cellpanel_idx", fields=['cell_panel']),
            models.Index(name="cell_source_idx", fields=['source']),
            # models.Index(name="org_taxid_idx", fields=['tax_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.cell_id} ({self.cell_line})"

    #------------------------------------------------
    @classmethod
    def str_CellID(cls,CellNo) -> str:
    #
    # Input:    CellClass GN, GP,...
    #           CellNo 
    # Output:   Cell_ID as string like GN_0001 
    #
        return(f"CL{CELL_SEP}{CellNo:04d}")

    #------------------------------------------------
    @classmethod
    def exists(cls,CellID,verbose=0):
    # Returns if an instance exists by cell_id
        return cls.objects.filter(cell_id=CellID).exists()

    #------------------------------------------------
    @classmethod
    def get(cls,CellID,verbose=0):
    # Returns an instance by cell_id
        try:
            retInstance = cls.objects.get(cell_id=CellID)
        except:
            if verbose:
                print(f"[CellID Not Found] {CellID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def find_Next_CellID(cls) -> str:
        Cell_IDSq = Sequence("Cell")
        Cell_nextID = next(Cell_IDSq)
        Cell_strID = cls.str_CellID(Cell_nextID)
        while cls.exists(Cell_strID):
            Cell_nextID = next(Cell_IDSq)
            Cell_strID = cls.str_CellID(Cell_nextID)
        return(Cell_strID)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.cell_id: 
            self.cell_id = self.find_Next_CellID()
            if self.cell_id: 
                super(Cell, self).save(*args, **kwargs)
        else:
            super(Cell, self).save(*args, **kwargs) 

#=================================================================================================
class Cell_Batch(AuditModel):
    """
    Cell Batch Collection
    """
#=================================================================================================

    HEADER_FIELDS = {
        "batch_id":"Batch ID",
        "batch_notes":"Batch Notes",
        "previous_batch_id":"Prev ID",
        "passage_number":"Passage",
        "qc_status":"QC",
        "batch_quality":"Batch Quality",
        "quality_source": "Quality by",
        "stock_date":"Stock Date",
        "stock_level":"Stock Levels",
        "biologist":"Biologist"
    }

    Choice_Dictionary = {
        'qc_status':'QC_Status',
        'batch_quality':'OrgBatch_Quality',
    }
    
    FORM_GROUPS = {
       'Group1': ["batch_id", "batch_notes", "qc_status", "batch_quality", "quality_source", "stock_date", "stock_level", "biologist" ]
       }
    #SEP = '_'

    alphanumeric = RegexValidator(r'^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')

    cellbatch_id  = models.CharField(primary_key=True, max_length=20, verbose_name = "CellBatch ID")
    cell_id = models.ForeignKey(Cell, null=False, blank=False, verbose_name = "Cell ID", on_delete=models.DO_NOTHING,
        db_column="cell_id", related_name="%(class)s_cell_id")
    previous_batch_id= models.CharField(max_length=20, verbose_name = "Previous CellBatch ID")
    passage_number= models.CharField(max_length=20, verbose_name = "Passage Number")
    batch_id  = models.CharField(max_length=12, null=False, blank=True, validators=[alphanumeric], verbose_name = "Batch ID")
    batch_notes= models.CharField(max_length=500, blank=True, verbose_name = "Batch Notes")
    batch_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Quality", on_delete=models.DO_NOTHING,
        db_column="batch_quality", related_name="%(class)s_batchquality")
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
        app_label = 'dcell'
        db_table = 'cellbatch'
        ordering=['cellbatch_id']
        indexes = [
            models.Index(name="cellbatch_cellbatch_idx",fields=['cell_id','batch_id']),
            models.Index(name="cellbatch_qc_idx",fields=['qc_status']),
            models.Index(name="cellbatch_sdate_idx",fields=['stock_date']),
            models.Index(name="cellbatch_slevel_idx",fields=['stock_level']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.cellbatch_id}"

    #------------------------------------------------
    @classmethod
    # Formats BatchNo:int -> BatchID:str 
    def str_BatchID(self,BatchNo:int) -> str:
        return(f"{BatchNo:02d}")
    #------------------------------------------------
    @classmethod
    # Formats CellID:str,BatchID:str -> CellBatchID:str
    def str_CellBatchID(self,CellID:str,BatchID:str) -> str:
        return(f"{CellID}{ORGBATCH_SEP}{BatchID}")

    #------------------------------------------------
    def find_Next_BatchID(self, CellID:str, BatchID:str=None) -> str:
        # Check for given BatchID    
        if BatchID:
            # Clean up BatchID - remove non alphanumeric character and make uppercase
            BatchID = re.sub(r'[^a-zA-Z0-9]', '', BatchID).upper()

            # Clean up BatchID - reformat numbers
            if BatchID.isnumeric():
                BatchID = self.str_BatchID(int(BatchID))

            next_CellBatch = self.str_CellBatchID(CellID,BatchID)
            if ~self.exists(next_CellBatch):
                return(BatchID)

        # Find new BatchID    
        next_BatchNo = 1
        next_CellBatch = self.str_CellBatchID(CellID,self.str_BatchID(next_BatchNo))
        while self.exists(next_CellBatch):
            next_BatchNo = next_BatchNo + 1
            next_CellBatch = self.str_CellBatchID(CellID,self.str_BatchID(next_BatchNo))
        return(self.str_BatchID(next_BatchNo))    

    #------------------------------------------------
    @classmethod
    def get(cls,CellBatchID,verbose=0):
    # Returns an instance if found by cellbatch_id
        try:
            retInstance = cls.objects.get(cellbatch_id=CellBatchID)
        except:
            if verbose:
                print(f"[CellBatch Not Found] {CellBatchID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,CellBatchID,verbose=0):
    # Returns if instance exists
        return cls.objects.filter(cellbatch_id=CellBatchID).exists()

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.cellbatch_id: 
            # creates new CellBatchID
            CellID = self.cell_id.cell_id
            BatchID = self.find_Next_BatchID(CellID,self.batch_id)
            if BatchID:
                self.batch_id = BatchID
                self.cellbatch_id = self.str_CellBatchID(CellID,BatchID)
                super(Cell_Batch,self).save(*args, **kwargs)
        else:
            # confirms Batch_ID from CellBatchID
            self.batch_id = str(self.cellbatch_id).replace(str(self.cell_id.cell_id),"").split(ORGBATCH_SEP)[1]
            super(Cell_Batch,self).save(*args, **kwargs)
        
# ================================================================================================


#=================================================================================================
class CellBatch_Stock(AuditModel):
    """
    Stock of Cell Batches
    
    """
#=================================================================================================
    HEADER_FIELDS={
        "cellbatch_id.cellbatch_id":{'CellBatch ID': {'cellbatch_id.cell_id.cell_id':LinkList["cell_id"]}},
        "cellbatch_id.cell_id.cell_name":"Cell",
        #"cellbatch_id.cell_id.cell_name":{'Cell ID': {'cell_id.cell_id.cell_id':LinkList['cell_id']}},
        "stock_type":"Stock Type",
        "n_created":"#Created",
        "n_left":"#Left",
        "location_freezer": "Freezer",
        "location_rack": "Rack",
        "location_column": "Column",
        "location_slot": "Slot",
        "stock_date":"Stock Date",
        "stock_note":"Stock Note",
        "biologist":"Biologist",
    }

    Choice_Dictionary = {
        'stock_type':'Stock_Type',
    }

    cellbatch_id = models.ForeignKey(Cell_Batch, null=False, blank=False, verbose_name = "CellBatch ID", on_delete=models.DO_NOTHING,
        db_column="cellbatch_id", related_name="%(class)s_cellbatch_id") 
    stock_type = models.ForeignKey(Dictionary, null=False, blank=False, verbose_name = "Stock Type", on_delete=models.DO_NOTHING,
        db_column="stock_type", related_name="%(class)s_stock")
    n_created = models.IntegerField(default=0, verbose_name = "#Created")
    n_left = models.IntegerField(default=0, verbose_name = "#Left")
    stock_date = models.DateField(null=True, blank = True, verbose_name = "Stock Date")
    stock_note = models.CharField(max_length=10, blank=True, verbose_name = "Stock Note")
    # passage_notes = models.CharField(max_length=30, blank=True, verbose_name = "Passage Notes")
    location_freezer = models.CharField(max_length=80, blank=True, verbose_name = "Freezer")
    location_rack = models.CharField(max_length=10, blank=True, verbose_name = "Rack")
    location_column = models.CharField(max_length=10, blank=True, verbose_name = "Column")
    location_slot = models.CharField(max_length=10, blank=True, verbose_name = "Slot")
    #stock_id = models.CharField(max_length=15, blank=True, verbose_name = "Stock ID")
    biologist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Biologist", on_delete=models.DO_NOTHING, 
        db_column="biologist", related_name="%(class)s_biologist")

    #------------------------------------------------
    class Meta:
        app_label = 'dcell'
        db_table = 'cellbatch_stock'
        ordering=['cellbatch_id','stock_type']
        indexes = [
            models.Index(name="cellbstock_stype_idx",fields=['stock_type']),
            models.Index(name="cellbstock_freezer_idx",fields=['location_freezer']),
            models.Index(name="cellbstock_stdate_idx",fields=['stock_date']),
            models.Index(name="cellbstock_nleft_idx",fields=['n_left']),

        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.pk} "
    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.cellbatch_id} {self.stock_type} {self.n_left}"

# ================================================================================================
