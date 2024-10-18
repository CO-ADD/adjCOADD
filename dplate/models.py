import re, math

from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from dscreen.models import Screen_Run
from dorganism.models import Organism_Batch
from dcell.models import Cell_Batch
from adjcoadd.constants import *



#-------------------------------------------------------------------------------------------------
# Plate Information - Test/Master Plates/Wells, and Labware
#-------------------------------------------------------------------------------------------------

PLATE_SIZE_DICT     = {24:(4,6), 48:(6,8), 96:(8,12), 384:(16,24), 1536:(32,48)}
PLATE_TYPE_DICT     = {'Plate':'Well','Rack':'Tube'}
PLATE_MATERIAL_DICT = {'PP':'Polypropylen','PS':'Polystyren','TC-PS':'TissueCulture','NBS-PS':'Non-Binding Surface'}

ROW_LABELS = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
                'Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF']

#=================================================================================================
class Labware(AuditModel):

    PLATE_SIZES = Choices(24,48,96,384,1536)
    PLATE_TYPES = Choices( ('Plate','Plate with Wells'),
                            ('Rack','Rack with Tubes')
                        )
    PLATE_MATERIAL = Choices('PP','PS','TC-PS','NBS-PS')
    PLATE_COLORS = Choices('Clear', 'Black', 'White')
    WELL_BOTTOMS = Choices('Clear', 'Black', 'White','Barcode')
    WELL_SHAPES = Choices('Flat','Round','U-Shape','V-Shape')
    
    
    Choice_Dictionary = {
        'plate_material':'Plate_Material',
    }

    labware_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Labware ID") 
    labware_name= models.CharField(max_length=50, blank=True, verbose_name = "Labware Name") 
    plate_size= models.PositiveSmallIntegerField(default=0, choices=PLATE_SIZES, blank=True, verbose_name = "Plate Size") 
    plate_type= models.CharField(max_length=10, choices=PLATE_TYPES, blank=True, verbose_name = "Plate Type") 
    plate_color = models.CharField(max_length=10, choices=PLATE_COLORS, blank=True, verbose_name = "Plate Type")
     
    #plate_material= models.CharField(max_length=8, choices=PLATE_MATERIAL, blank=True, verbose_name = "Material") 
    plate_material= models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Material", on_delete=models.DO_NOTHING,
        db_column="material", related_name="%(class)s_material")
    well_shape= models.CharField(max_length=10, choices=WELL_SHAPES, verbose_name = "Shape") 
    well_bottom= models.CharField(max_length=10, choices=WELL_BOTTOMS, blank=True, verbose_name = "Bottom") 
    working_volume = models.DecimalField(max_digits=9, decimal_places=1, default=0, blank=True, verbose_name = "Working volume (uL)")
    
    brand= models.CharField(max_length=25, blank=True, verbose_name = "Brand") 
    model= models.CharField(max_length=25, blank=True, verbose_name = "Model") 

    class Meta:
        app_label = 'dplate'
        db_table = 'labware'
        ordering=['labware_id']

#=================================================================================================
class Well(AuditModel):
    """
    An abstract Well class model that provides general Well properties/method 
    """
    PLATE_CLASS = Plate

    #------------------------------------------------
    plate_id = models.ForeignKey(PLATE_CLASS,  verbose_name = "Plate ID", on_delete=models.DO_NOTHING,
        db_column="plate_id", related_name="%(class)s_plateid")
    well_id = models.CharField(max_length=6, verbose_name = "Well ID")

    class Meta:
        abstract = True
        ordering=['plate_id','well_id']
        indexes = [
            models.Index(fields=['plate_id']),
            models.Index(fields=['well_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.plate_id} {self.well_id}"
    #------------------------------------------------
    def __repr__(self) -> str:
        # return f"{self.__name__}: {self.pk}"
        return f"{self.plate_id} {self.well_id}"

    #------------------------------------------------
    @classmethod
    def get(cls,PlateID, WellID, verbose=0):
        try:
            retInstance = cls.objects.get(plate_id=PlateID, well_id=WellID)
        except:
            if verbose:
                print(f"[Well Not Found] {PlateID} {WellID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,PlateID,WellID,verbose=0):
        return cls.objects.filter(plate_id=PlateID, well_id=WellID).exists()


#=================================================================================================


#=================================================================================================
class Plate(AuditModel):
    """
    An abstract Plate class model that provides general Plate properties/method 
    """
#=================================================================================================

    WELL_CLASS = Well
    
    PLATE_SIZES = {24:(4,6), 48:(6,8), 96:(8,12), 384:(16,24), 1536:(32,48)}
    ROW_LABELS = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
                 'Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF']
    MAP_POSITIONS = {'wellID':0,'pos2D':1,'pos1D':2}

    Choice_Dictionary = {
        'plate_type':'Plate_Type',
    }

    plate_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Plate ID")
    plate_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Type", on_delete=models.DO_NOTHING,
        db_column="plate_type", related_name="%(class)s_platetype")
    labware_id = models.ForeignKey(Labware, null=True, blank=True, verbose_name = "Labware ID", on_delete=models.DO_NOTHING,
        db_column="labware_id", related_name="%(class)s_labwareid")
    rows =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "Rows") 
    cols =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "Cols") 
    n_wells =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "nWells")

    wells = {}

    #------------------------------------------------
    class Meta:
        abstract = True
        ordering=['plate_id']
        indexes = [
            models.Index(fields=['plate_id']),
        ]


    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.plate_id}"
    #------------------------------------------------
    def __repr__(self) -> str:
        # return f"{self.__name__}: {self.pk}"
        return f"{self.plate_id}"

    #------------------------------------------------
    @classmethod
    def get(cls,PlateID,WellData=True,verbose=0):
        try:
            retInstance = cls.objects.get(plate_id=PlateID)
            if WellData:
                retInstance.n_wells = retInstance.get_wells()
        except:
            if verbose:
                print(f"[Plate Not Found] {PlateID} ")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,PlateID,verbose=0):
        return cls.objects.filter(plate_id=PlateID).exists()

    #------------------------------------------------
    def set_platesize(self,PlateSize):
        # -- Set Plate Size
        if isinstance(PlateSize,int):
            if PlateSize in self.PLATE_SIZES:
                self.rows = self.PLATE_SIZES[PlateSize][0]
                self.cols = self.PLATE_SIZES[PlateSize][1]
            else:
                raise KeyError(f"Undefined PlateSize {PlateSize}")
        elif isinstance(PlateSize,tuple):
            self.rows = PlateSize[0]
            self.cols = PlateSize[1]
        else:
            self.rows = 0
            self.cols = 0
            raise KeyError(f"Undefined PlateSize parameters {PlateSize}")
        self.plate_size = self.rows * self.cols

    #------------------------------------------------
    @classmethod
    def create(cls,PlateID,PlateSize,PlateType,WellData=True):
        _plate = cls()
        _plate.set_platesize(PlateSize)
        _plate.plate_type = Dictionary.get(cls.Choice_Dictionary["plate_type"],PlateType)
        _plate.create_wells()
        return(_plate)
    
    #------------------------------------------------
    def get_wells(self) -> int:
        pass

    #------------------------------------------------
    def save_wells(self) :
        pass

    #------------------------------------------------
    def map_well(self,loc,check=True):
        if isinstance(loc,int) :
            if check:
                if not loc in self.validate['pos1D']:
                    raise Exception(f"Invalid 1D Plate Position: {loc}")
            row = int(math.ceil(float(loc)/float(self.cols))) - 1
            col = loc - (row * self.cols) - 1
        elif isinstance(loc,tuple):
            if check:
                if not loc in self.validate['pos2D']:
                    raise Exception(f"Invalid 2D Plate Position: {loc}")
            row = loc[0] - 1
            col = loc[1] - 1
        elif isinstance(loc,str) :
            res = re.findall('([A-Za-z]+|\d+)',loc)
            loc = f"{res[0]}{int(res[1]):02d}"
            if check:
                if not loc in self.validate['wellID']:
                    raise Exception(f"Invalid Well ID: {loc}")
            row = self.ROW_LABELS.index(loc[0])
            col = int(loc[1:]) - 1
        else:
            raise  Exception(f"Unrecognized Plate Location Type: {loc}")

        pos = self.cols * row + col +1
        id = f"{self.ROW_LABELS[row]:s}{(col+1):02d}"
        return(id,(row+1,col+1),pos)

    def well_pos(self,loc):
        m = self.map_wellmap(loc)
        return(m[3])

    def well_rowcol(self,loc):
        m = self.map_well(loc)
        return(m[1])

    def well_id(self,loc):
        m = self.map_well(loc)
        return(m[0])

#=================================================================================================
class TestPlate(Plate):
    """

    """
#=================================================================================================

    WELL_CLASS = TestWell
    RESULT_TYPES = Choices('MIC','CC50','HC50','SYN-MIC')

    Choice_Dictionary = {
        'result_type':'Result_Type',
        'plate_quality':'Plate_Quality',
    }

    # plate_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Plate ID")
    # labware_id = models.ForeignKey(Labware, null=True, blank=True, verbose_name = "Labware ID", on_delete=models.DO_NOTHING,
    #     db_column="labware_id", related_name="%(class)s_labwareid")
    
    motherplate_ids = ArrayField(models.CharField(max_length=25, null=True, blank=True), size=2, verbose_name = "Mother Plates", null=True, blank=True)
                                
    # Plate_Size = models.CharField(max_length=10)
    # N_Wells = models.PositiveIntegerField()
    prep_date = models.DateField(null=True, blank=True, verbose_name = "Prep Date")
    plating = models.CharField(max_length=10, blank=True, verbose_name = "Plating by")
    run_id = models.ForeignKey(Screen_Run, null=True, blank=True, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_runid")    
    result_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Result Type", on_delete=models.DO_NOTHING,
        db_column="result_type", related_name="%(class)s_resulttype")
    
    assay_id = models.CharField(max_length=25, blank=True, verbose_name = "Assay ID")
    test_date = models.DateField(null=True, blank=True, verbose_name = "Test Date")
    test_media = models.CharField(max_length=50, blank=True, verbose_name = "Media")
    #test_strain = models.CharField(max_length=15, blank=True, verbose_name = "Strain")

    test_strain = models.ForeignKey(Organism_Batch, null=True, blank=True, verbose_name = "Strain", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatchid")
    test_cell = models.ForeignKey(Cell_Batch, null=True, blank=True, verbose_name = "Cell", on_delete=models.DO_NOTHING,
        db_column="cellbatch_id", related_name="%(class)s_cellbatchid")
    
    test_dye = models.CharField(max_length=25, blank=True, verbose_name = "Dye")
    test_addition = models.CharField(max_length=25, blank=True, verbose_name = "Addition")
    test_volume = models.DecimalField(max_digits=10, decimal_places=2, verbose_name = "Volume (uL)")
    test_processing = models.CharField(max_length=25, blank=True, verbose_name = "Processing")
    test_issues = models.CharField(max_length=100, blank=True, verbose_name = "Issue")

    reader = models.CharField(max_length=50, blank=True, verbose_name = "Reader")
    n_reads = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "#Reads")
    readout_type = models.CharField(max_length=25, blank=True, verbose_name = "Readout Type")
    experiment = models.CharField(max_length=80, blank=True, verbose_name = "Experiment")
    protocol = models.CharField(max_length=80, blank=True, verbose_name = "Protocol")
    input_file = models.CharField(max_length=80, blank=True, verbose_name = "Input File")
    test_operator = models.CharField(max_length=100, blank=True, verbose_name = "Operator")
    
    has_readout = models.BooleanField(default=False, verbose_name = "Has Readout")
    has_sample = models.BooleanField(default=False, verbose_name = "Has Sample")
    has_layout = models.BooleanField(default=False, verbose_name = "Has Layout")
    has_inhibition = models.BooleanField(default=False, verbose_name = "Has Inhibition")
    has_doseresponse = models.BooleanField(default=False, verbose_name = "Has Doseresponse")
    process_status = models.PositiveSmallIntegerField(default=0, verbose_name = "Process Status")

    control_layout = models.CharField(max_length=25, blank=True, verbose_name = "Layout")
    synergy_samples =ArrayField(models.CharField(max_length=100, null=True, blank=True), size=2, verbose_name = "Synergy Samples", null=True, blank=True)

    # Control_Count = models.PositiveIntegerField()
    # Layout_Dilution = models.CharField(max_length=25)

    positive_control = ArrayField(models.DecimalField(max_digits=7, decimal_places=2),size=4)
    negative_control = ArrayField(models.DecimalField(max_digits=7, decimal_places=2),size=4)
    sample_stats = ArrayField(models.DecimalField(max_digits=7, decimal_places=2),size=4)
    edge_stats = ArrayField(models.DecimalField(max_digits=7, decimal_places=2),size=2)
    analysis_parameter = models.CharField(max_length=100, blank=True, verbose_name = "Analysis")
    zfactor = models.DecimalField(max_digits=7, decimal_places=2)
    plate_qc = models.DecimalField(max_digits=7, decimal_places=2)
    plate_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Quality", on_delete=models.DO_NOTHING,
        db_column="plate_quality", related_name="%(class)s_platequality")
    
    # Test_Dye_Conc = models.DecimalField(max_digits=10, decimal_places=2)
    # Test_Dye_Conc_Unit = models.CharField(max_length=10)
    # Signal_Window = models.DecimalField(max_digits=7, decimal_places=2)
    # NegControl_Median = models.DecimalField(max_digits=10, decimal_places=2)
    # NegControl_MAD = models.DecimalField(max_digits=10, decimal_places=2)
    # PosControl_Median = models.DecimalField(max_digits=10, decimal_places=2)
    # PosControl_MAD = models.DecimalField(max_digits=10, decimal_places=2)
    # Sample_Median = models.DecimalField(max_digits=10, decimal_places=2)
    # Sample_MAD = models.DecimalField(max_digits=10, decimal_places=2)
    # Edge_Median = models.DecimalField(max_digits=7, decimal_places=2)
    # NonEdge_Median = models.DecimalField(max_digits=7, decimal_places=2)

    class Meta:
        app_label = 'dplate'
        db_table = 'testplate'
        ordering=['run_id','plate_id']
        indexes = [
            models.Index(name="testplate_labw_idx",fields=['labware_id']),
            # models.Index(fields=['MotherPlate_ID']),
            models.Index(name="testplate_rest_idx",fields=['result_type']),
            models.Index(name="testplate_ass_idx",fields=['assay_id']),
            models.Index(name="testplate_run_idx",fields=['run_id']),
            models.Index(name="testplate_read_idx",fields=['readout_type']),
            models.Index(name="testplate_proc_idx",fields=['process_status']),
            models.Index(name="testplate_qc_idx",fields=['plate_qc']),
            models.Index(name="testplate_pq_idx",fields=['plate_quality']),
            models.Index(name="testplate_has_idx",fields=['has_readout', 'has_sample', 'has_layout', 'has_inhibition', 'has_doseresponse']),
            models.Index(name="testplate_test_idx",fields=['test_media', 'test_strain', 'test_dye', 'test_addition']),
        ]
        
#=================================================================================================
class TestWell(AuditModel):
    """

    """
#=================================================================================================

    RESULT_TYPES = Choices('MIC','CC50','HC50','SYN-MIC')

    Choice_Dictionary = {
        'solvent_conc_unit':'Conc_Unit'
    }

    plate_id = models.ForeignKey(TestPlate, null=True, blank=True, verbose_name = "Plate ID", on_delete=models.DO_NOTHING,
        db_column="plate_id", related_name="%(class)s_plateid")
    well_id = models.CharField(max_length=5, blank=False, verbose_name = "Well ID")

    n_samples = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "#Samples")    
    samplebatch_id = ArrayField(models.CharField(max_length=25),size=4)
    analysis_set = ArrayField(models.CharField(max_length=5),size=4)
    sample_concs = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),size=4)
    conc_units = ArrayField(models.CharField(max_length=5),size=4)
    conc_types = ArrayField(models.CharField(max_length=5),size=4)
    solvent = models.CharField(max_length=25)
    solvent_conc = models.DecimalField(max_digits=12, decimal_places=4)
    solvent_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Conc Unit", on_delete=models.DO_NOTHING,
         db_column="solvent_conc_unit", related_name="%(class)s_solvconcunit")

    is_control = models.BooleanField(default=False, verbose_name = "is Control")
    is_poscontrol = models.BooleanField(default=False, verbose_name = "is PosCtrl")
    is_negcontrol = models.BooleanField(default=False, verbose_name = "is NegCtrl")
    is_sample = models.BooleanField(default=False, verbose_name = "is Sample")
    is_skip = models.BooleanField(default=False, verbose_name = "is Skip")
    is_valid = models.BooleanField(default=False, verbose_name = "is Valid")
    
    readouts = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),size=4)
    readout_types = ArrayField(models.CharField(max_length=5),size=4)
    
    inhibition = models.DecimalField(max_digits=7, decimal_places=2)
    zscore = models.DecimalField(max_digits=7, decimal_places=2)
    mscore = models.DecimalField(max_digits=7, decimal_places=2)
    bscore = models.DecimalField(max_digits=7, decimal_places=2)
    
    class Meta:
        app_label = 'dplate'
        db_table = 'testwell'
        ordering=['plate_id','well_id']
        indexes = [
            models.Index(name="testwell_tptw_idx",fields=['plate_id','well_id']),
            models.Index(name="testwell_smpid_idx",fields=['samplebatch_id']),
            models.Index(name="testwell_nsmp_idx",fields=['n_samples']),
            models.Index(name="testwell_is_idx",fields=['is_control', 'is_poscontrol', 'is_negcontrol', 'is_sample']),
            models.Index(name="testwell_skip_idx",fields=['is_skip', 'is_valid']),
        ]
