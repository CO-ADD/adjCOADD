import re, math
import pandas as pd
import numpy as np

from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify
from django.forms.models import model_to_dict

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from apputil.utils.data import addto_StrList, strList_to_List
from dscreen.models import Screen_Run, Assay
from dorganism.models import Organism_Batch
from dcell.models import Cell_Batch
from dsample.models import Sample_Base
from ddrug.utils.bio_data import ActScoreSC_Cutoff, ActType_SC
from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

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
    #PLATE_MATERIAL = Choices('PP','PS','TC-PS','NBS-PS')
    PLATE_COLORS = Choices('Clear', 'Black', 'White')
    WELL_BOTTOMS = Choices('Clear', 'Black', 'White','Barcode')
    WELL_SHAPES = Choices('Flat','Round','U-Shape','V-Shape')
    WELL_TYPE = Choices('Well','Tube')
    WELL_SIZE = Choices('Shallow','Deep','Storage')
    
    
    Choice_Dictionary = {
        'plate_material':'Plate_Material',
    }

    labware_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Labware ID") 
    labware_name= models.CharField(max_length=50, blank=True, verbose_name = "Labware Name") 
    labware_type= models.CharField(max_length=15, blank=True, verbose_name = "Labware Type") 
    labware_notes= models.CharField(max_length=50, blank=True, verbose_name = "Labware Type") 
    plate_size= models.PositiveSmallIntegerField(default=0, choices=PLATE_SIZES, blank=True, verbose_name = "Plate Size") 
    plate_type= models.CharField(max_length=10, choices=PLATE_TYPES, blank=True, verbose_name = "Plate Type") 
    plate_color = models.CharField(max_length=10, choices=PLATE_COLORS, blank=True, verbose_name = "Plate Type")
     
    plate_material= models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Material", on_delete=models.DO_NOTHING,
        db_column="material", related_name="%(class)s_material")
    well_shape= models.CharField(max_length=20, choices=WELL_SHAPES, verbose_name = "Shape") 
    well_bottom= models.CharField(max_length=20, choices=WELL_BOTTOMS, blank=True, verbose_name = "Bottom") 
    well_type= models.CharField(max_length=20, choices=WELL_TYPE, verbose_name = "Type") 
    well_size= models.CharField(max_length=20, choices=WELL_SIZE, blank=True, verbose_name = "Size") 
    working_volume = models.DecimalField(max_digits=9, decimal_places=1, default=0, blank=True, verbose_name = "Working volume (uL)")
    
    brand= models.CharField(max_length=25, blank=True, verbose_name = "Brand") 
    model= models.CharField(max_length=25, blank=True, verbose_name = "Model") 

    class Meta:
        app_label = 'dplate'
        db_table = 'labware'
        ordering=['labware_id']

#=================================================================================================



#=================================================================================================
class Plate(AuditModel):
    """
    An abstract Plate class model that provides general Plate properties/method 
    """
#=================================================================================================

    #WELL_CLASS = Well
    
    PLATE_SIZES = {24:(4,6), 48:(6,8), 96:(8,12), 384:(16,24), 1536:(32,48)}
    ROW_LABELS = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
                 'Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF']
    MAP_POSITIONS = {'wellID':0,'pos2D':1,'pos1D':2}
    WELL_POS = MAP_POSITIONS['wellID']

    Choice_Dictionary = {
        'plate_type':'Plate_Type',
    }

    plate_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Plate ID")
    plate_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Type", on_delete=models.DO_NOTHING,
        db_column="plate_type", related_name="%(class)s_platetype")
    labware_id = models.ForeignKey(Labware, null=True, blank=True, verbose_name = "Labware ID", on_delete=models.DO_NOTHING,
        db_column="labware_id", related_name="%(class)s_labwareid")
    n_rows =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "nRows") 
    n_cols =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "nCols") 
    n_wells =  models.PositiveSmallIntegerField(default=0, blank=True, verbose_name = "nWells")

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
        except:
            if verbose:
                logger.warning(f"[Plate Not Found] {PlateID} ")
            retInstance = None

        if retInstance and WellData:
            retInstance.n_wells = retInstance.get_wells()

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
                self.n_rows = self.PLATE_SIZES[PlateSize][0]
                self.n_cols = self.PLATE_SIZES[PlateSize][1]
                self.n_wells = PlateSize
            else:
                raise KeyError(f"Undefined PlateSize {PlateSize}")
        elif isinstance(PlateSize,tuple):
            self.n_rows = PlateSize[0]
            self.n_cols = PlateSize[1]
            self.n_wells = self.n_rows * self.n_cols
        else:
            self.n_rows = 0
            self.n_cols = 0
            self.n_wells = 0
            raise KeyError(f"Undefined PlateSize parameters {PlateSize}")

    #------------------------------------------------
    @classmethod
    def new(cls,PlateID,PlateSize,PlateType,WellData=True):
        #print(f"[Plate.new] {PlateID} {PlateSize} {PlateType}")
        _plate = cls()
        _plate.plate_id = PlateID.upper()
        _plate.set_platesize(PlateSize)
        _plate.plate_type = Dictionary.get(cls.Choice_Dictionary["plate_type"],PlateType)
        _plate.init_wells()
        return(_plate)
    
    #------------------------------------------------
    # Well Mapping
    #------------------------------------------------
    def map_well(self,loc,check=True):
        if isinstance(loc,int) :
            if check:
                if not loc in self.well_check['pos1D']:
                    raise Exception(f"{self.plate_id} Invalid 1D Plate Position: {loc}")
            row = int(math.ceil(float(loc)/float(self.n_cols))) - 1
            col = loc - (row * self.n_cols) - 1
        elif isinstance(loc,tuple):
            if check:
                if not loc in self.well_check['pos2D']:
                    raise Exception(f"{self.plate_id} Invalid 2D Plate Position: {loc}")
            row = loc[0] - 1
            col = loc[1] - 1
        elif isinstance(loc,str) :
            res = re.findall('([A-Za-z]+|\d+)',loc)
            loc = f"{res[0]}{int(res[1]):02d}"
            if check:
                if not loc in self.well_check['wellID']:
                    raise Exception(f"{self.plate_id} Invalid Well ID: {loc}")
            row = self.ROW_LABELS.index(loc[0])
            col = int(loc[1:]) - 1
        else:
            raise  Exception(f"Unrecognized Plate Location Type: {loc}")

        pos = self.n_cols * row + col +1
        id = f"{self.ROW_LABELS[row]:s}{(col+1):02d}"
        return(id,(row+1,col+1),pos)

    #------------------------------------------------
    def well_pos(self,loc):
        m = self.map_well(loc)
        return(m[self.MAP_POSITIONS['pos1D']])

    #------------------------------------------------
    def well_rowcol(self,loc):
        m = self.map_well(loc)
        return(m[self.MAP_POSITIONS['pos2D']])

    #------------------------------------------------
    def well_id(self,loc):
        m = self.map_well(loc)
        return(m[self.MAP_POSITIONS['wellID']])

    #------------------------------------------------
    # Initialising Wells
    #------------------------------------------------
    def get_wells(self) -> int:
        # Create empty Wells
        self.init_wells()
        # Fill with Database Wells
        # to be implmented in specific models

    #------------------------------------------------
    def save_wells(self) :
        if self.wells:
            for w in self.wells:
                if self.wells[w] is not None:
                    self.wells[w].save()
        
    #------------------------------------------------
    def is_edgewell(self,well_id):
        (r,c) = self.well_rowcol(well_id)
        return(r == 1 or r == self.n_rows or c == 1 or c == self.n_cols)

    #--------------------------------------------------------------
    def init_wells(self, WellModel=None, PlateInstance = None, reset=False) -> int:
    #
    # Initalise Wells Dictionary, with Empty or with New WellModels() 
    #
        #print(f"[Plate.init_wells] {WellModel} {reset}")
        if not hasattr(self,'wells') or reset:

            # Well Check Arrays
            self.well_check = {}
            self.wells = {}

            for key in self.MAP_POSITIONS:
                self.well_check[key] = []

            # Create Dict of Wells and Well_Check
            for n in range(1,self.n_wells+1):
                m = self.map_well(n,check=False)
                if WellModel is not None and PlateInstance is not None:
                    #print(f"[Plate.init_wells] {m} with {WellModel} for {PlateInstance}")

                    self.wells[m[self.WELL_POS]] = WellModel()
                    self.wells[m[self.WELL_POS]].well_id = m[0]
                    self.wells[m[self.WELL_POS]].plate_id = PlateInstance
                else:
                    self.wells[m[self.WELL_POS]] = None

                for key in self.MAP_POSITIONS:
                    self.well_check[key].append(m[self.MAP_POSITIONS[key]])
    
    #------------------------------------------------
    def fill_wells(self, WellModel=None):
        if self.wells:
            for w in self.wells:
                if self.wells[w] is None:
                    self.wells[w] = WellModel()
                    self.wells[w].well_id = w
                    self.wells[w].plate_id = self.plate_id

    #------------------------------------------------
    def reset_well_fields(self, Fields=[]):
        if Fields:
            if self.wells:
                for w in self.wells:
                    if self.wells[w]:
                        for f in Fields:
                            self.wells[w].init_field(f, Reset=True)

    #--------------------------------------------------------------
    def get_well_fielddata(self, Field, Selection = None):
    #
    # Get fielddata from the Wells, by Selection
    # Model specific implemnetation    
        pass

#=================================================================================================#=================================================================================================
class TestPlate(Plate):
    """

    """
#=================================================================================================

    #WELL_CLASS = TestWell
    RESULT_TYPES = Choices('MIC','CC50','HC50','SYN-MIC')
    ZFACTOR_CUTOFF = 0.2

    STATS_MEDIAN = 0
    STATS_MAD = 1
    STATS_MEAN = 2
    STATS_STD = 3

    Choice_Dictionary = {
        'result_type':'Result_Type',
        'plate_quality':'Data_Quality',
        'plate_type':'Plate_Type',
    }

    # plate_id = models.CharField(primary_key=True, max_length=25, verbose_name = "Plate ID")
    # labware_id = models.ForeignKey(Labware, null=True, blank=True, verbose_name = "Labware ID", on_delete=models.DO_NOTHING,
    #     db_column="labware_id", related_name="%(class)s_labwareid")
    
    motherplate_ids = ArrayField(models.CharField(max_length=25, null=True, blank=True), size=2, verbose_name = "Mother Plates", null=True, blank=True)
                                
    # Plate_Size = models.CharField(max_length=10)
    # N_Wells = models.PositiveIntegerField()
    # prep_date = models.DateField(null=True, blank=True, verbose_name = "Prep Date")

    plating = models.CharField(max_length=10, blank=True, verbose_name = "Plating by")
    run_id = models.ForeignKey(Screen_Run, null=True, blank=True, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_runid")    
    result_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Result Type", on_delete=models.DO_NOTHING,
        db_column="result_type", related_name="%(class)s_resulttype")
    
    assay_id = models.ForeignKey(Assay, null=True, blank=True, verbose_name = "Assay ID", on_delete=models.DO_NOTHING,
        db_column="assay_id", related_name="%(class)s_assay_id")
        
    ora_assay_id = models.CharField(max_length=25, blank=True, verbose_name = "ora Assay ID")
    test_date = models.DateField(null=True, blank=True, verbose_name = "Test Date")
    #test_strain = models.CharField(max_length=15, blank=True, verbose_name = "Strain")

    test_orgbatch = models.ForeignKey(Organism_Batch, null=True, blank=True, verbose_name = "OrgBatch", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatchid")
    test_cellbatch = models.ForeignKey(Cell_Batch, null=True, blank=True, verbose_name = "CellBatch", on_delete=models.DO_NOTHING,
        db_column="cellbatch_id", related_name="%(class)s_cellbatchid")
    
    # Specific Assay condition - extension of Assay Condition
    test_media = models.CharField(max_length=50, blank=True, verbose_name = "Media")
    test_dye = models.CharField(max_length=25, blank=True, verbose_name = "Dye")
    test_additive = models.CharField(max_length=25, blank=True, verbose_name = "Additive")
    subculture_type = models.CharField(max_length=25, blank=True, verbose_name = "Subculture/Seeding")
    incubation_time = models.CharField(max_length=25, blank=True, verbose_name = "Incubation Time")

    test_volume = models.DecimalField(max_digits=10, decimal_places=2, default = -1, verbose_name = "Volume (uL)")
    test_processing = models.CharField(max_length=25, blank=True, verbose_name = "Processing")
    test_issues = models.CharField(max_length=150, blank=True, verbose_name = "Issue")

    reader = models.CharField(max_length=50, blank=True, verbose_name = "Reader")
    n_readouts = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Readouts")
    readout_type = models.CharField(max_length=25, blank=True, verbose_name = "Readout Type")
    experiment = models.CharField(max_length=80, blank=True, verbose_name = "Experiment")
    protocol = models.CharField(max_length=80, blank=True, verbose_name = "Protocol")
    input_file = models.CharField(max_length=80, blank=True, verbose_name = "Input File")
    test_operator = models.CharField(max_length=100, blank=True, verbose_name = "Operator")
    
    n_reads = models.SmallIntegerField(default=-1,verbose_name = "#Reads")
    n_samples = models.SmallIntegerField(default=-1, verbose_name = "#Samples")
    n_layout = models.SmallIntegerField(default=-1, verbose_name = "Layout")
    n_inhibitions = models.SmallIntegerField(default=-1, verbose_name = "#Inhibition")
    n_doseresponses = models.SmallIntegerField(default=-1, verbose_name = "#Doseresponse")
    n_synergies = models.SmallIntegerField(default=-1, verbose_name = "#Synergies")
    process_status = models.SmallIntegerField(default=-1, verbose_name = "Process Status")

    control_layout = models.CharField(max_length=35, blank=True, verbose_name = "Layout")
    synergy_cmpbatches =ArrayField(models.CharField(max_length=100, null=True, blank=True), 
                                   size=2, verbose_name = "Synergy CmpBatches", null=True, blank=True)

    # Control_Count = models.PositiveIntegerField()
    # Layout_Dilution = models.CharField(max_length=25)

    poscontrol_stats = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),
                                  size=4, null=True, blank=True, verbose_name = "PosCtrl")
    negcontrol_stats = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),
                                  size=4, null=True, blank=True, verbose_name = "NegCtrl")
    sample_stats = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),
                              size=4, null=True, blank=True, verbose_name = "Sample")
    edge_stats = ArrayField(models.DecimalField(max_digits=12, decimal_places=4),
                            size=2, null=True, blank=True, verbose_name = "Edge")
    
    analysis = models.CharField(max_length=100, blank=True, verbose_name = "Analysis")
    zfactor = models.DecimalField(max_digits=12, decimal_places=3, default=-1,)
    plate_qc = models.DecimalField(max_digits=12, decimal_places=3, default=-1)
    plate_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Quality", on_delete=models.DO_NOTHING,
        db_column="plate_quality", related_name="%(class)s_platequality")
    
    # Test_Dye_Conc = models.DecimalField(max_digits=10, decimal_places=2)
    # Test_Dye_Conc_Unit = models.CharField(max_length=10)
    # Signal_Window = models.DecimalField(max_digits=7, decimal_places=2)

    class Meta:
        app_label = 'dplate'
        db_table = 'testplate'
        ordering=['run_id','plate_id']
        indexes = [
            models.Index(name="testplate_labw_idx",fields=['labware_id']),
            # models.Index(fields=['MotherPlate_ID']),
            models.Index(name="testplate_rest_idx",fields=['result_type']),
            models.Index(name="testplate_oass_idx",fields=['ora_assay_id']),
            models.Index(name="testplate_run_idx",fields=['run_id']),
            models.Index(name="testplate_read_idx",fields=['readout_type']),
            models.Index(name="testplate_proc_idx",fields=['process_status']),
            models.Index(name="testplate_qc_idx",fields=['plate_qc']),
            models.Index(name="testplate_pq_idx",fields=['plate_quality']),
            models.Index(name="testplate_nnn_idx",fields=['n_reads', 'n_samples', 'n_layout', 'n_inhibitions', 'n_doseresponses','n_synergies']),
        #    models.Index(name="testplate_test_idx",fields=['test_media', 'test_strain', 'test_dye', 'test_addition']),
        ]

    #--------------------------------------------------------------
    def __repr__(self):
        _str  = f" [Testplate] {self.plate_id} Size:{self.n_wells} "
        _str += f"[R:{self.n_reads} S:{self.n_samples} L:{self.n_layout}"
        _str += f"I:{self.n_inhibitions} DR:{self.n_doseresponses} SYN:{self.n_synergies}]"
        if hasattr(self,'wells'):
            _wellid = list(self.wells.keys())
            _str += f" Wells:{len(self.wells)} [{_wellid[0]}..{_wellid[-1]}]"
           
        return(_str)
    
    #------------------------------------------------
    @classmethod
    def new(cls,PlateID,PlateSize,WellData=True):
        #print(f"[TestPlate.new] {PlateID} {PlateSize} ")
        if cls.exists(PlateID.upper()):
            logger.warning(f"[Testplate] New {PlateID.upper()} alreday exists ")
            return(None)
        else:
            _plate = cls()
            _plate.plate_id = PlateID.upper()
            _plate.set_platesize(PlateSize)
            _plate.plate_type = Dictionary.get(cls.Choice_Dictionary["plate_type"],'Test')

            if WellData:
                #print(f"[TestPlate.new] WithModel {TestWell}")
                _plate.init_wells(WellModel=TestWell,PlateInstance=_plate)
            else:
                _plate.init_wells(WellModel=None, PlateInstance=None)
            return(_plate)

    #------------------------------------------------
    def validate_model(self, WellData=True, verbose = 0):
        retDict = []
        PlateDict = super(TestPlate, self).validate_model(verbose=verbose)
        for wd in PlateDict:
            retDict.append(wd)

        if self.wells and WellData:
            for w in self.wells:
                if self.wells[w] is not None:
                    WellDict = super(TestWell,self.wells[w]).validate_model(verbose=verbose)
                    for wd in WellDict:
                        if 'plate_id' not in wd :
                            retDict.append(wd)
        
        return(retDict)

    #------------------------------------------------
    def init_model(self, WellData=True, verbose = 0):
        retDict = []
        super(TestPlate, self).init_fields()

        if self.wells and WellData:
            for w in self.wells:
                if self.wells[w] is not None:
                    super(TestWell,self.wells[w]).init_fields()
        
    #------------------------------------------------
    def save(self, *args, **kwargs):
        if self.plate_id:
            verbose = kwargs.get('verbose',0)
            kwargs.pop("verbose",None)

            super(TestPlate, self).save(*args, **kwargs)
            if verbose > 0:
                logger.info(f"[TestPlate.save] {self.plate_id}")
            if hasattr(self,'wells'):
                for w in self.wells:
                    if self.wells[w] is not None:
                        if verbose > 1:
                            logger.info(f"[TestWell.save] {self.wells[w]}")        
                        super(TestWell,self.wells[w]).save(*args, **kwargs)
        else:
            logger.warning(f"[TestPlate] SAVE has no PlateID ")

    #--------------------------------------------------------------
    def set_well_field(self,pos,field,value):
        if hasattr(self,'wells'):
            _w = self.well_id(pos)
            setattr(self.wells[_w],field,value)

    #--------------------------------------------------------------
    def get_wells(self, fill_missing=True) -> int:
        # Create None Wells
        self.init_wells(WellModel=None, PlateInstance=None)

        # get Wells for that Plate 
        qryTW = TestWell.objects.filter(plate_id=self)
        lWells = qryTW.count()
        for w in qryTW:
            m = self.map_well(w.well_id)
            self.wells[m[0]] = w

        # Fill None Wells with empty TestWell
        if lWells < len(self.wells) and fill_missing:
            for w in self.wells:
                if self.wells[w] is None:
                    self.wells[w] = TestWell()
                    self.wells[w].well_id = w
                    self.wells[w].plate_id = self
            lWells = len(self.wells)
        
        return(lWells)

    #--------------------------------------------------------------
    def get_welldata(self,) -> pd.DataFrame:
        _dicts = []
        if self.wells:
            for w in self.wells:
                if self.wells[w] is not None:        
                    _dicts.append(self.wells[w].get_welldict())
        self.well_data = pd.DataFrame(_dicts)
        return(len(self.well_data))

    #--------------------------------------------------------------
    def get_readouts(self,Selection = None):
        _readouts = []
        if self.wells:
            for w in self.wells:
                if self.wells[w] is not None:
                    if Selection:
                        if Selection == 'is_edge':
                            if self.is_edgewell(w) and not getattr(self.wells[w],'is_negcontrol'):
                                _readouts.append(float(self.wells[w].readouts[0]))
                        elif Selection == 'is_nonedge':
                            if not self.is_edgewell(w) and not getattr(self.wells[w],'is_negcontrol'):
                                _readouts.append(float(self.wells[w].readouts[0]))
                        else:
                            if getattr(self.wells[w],Selection):
                                _readouts.append(float(self.wells[w].readouts[0]))
                    else:
                        _readouts.append(float(self.wells[w].readouts[0]))
        return(np.array(_readouts))



    #--------------------------------------------------------------
    def get_well_fielddata(self, Field, Selection = None):
    #
    # Get fielddata from the Wells, by Selection
    # Model specific implemnetation    
        pass

    #--------------------------------------------------------------
    def apply_layout(self,verbose=0) -> int:
        CONTROL_LABELS = ['is_negcontrol','is_poscontrol','is_control','is_sample']
        CONTROL_ORDER = {'Neg':['is_negcontrol'],'Pos':['is_poscontrol'],'Ref':['is_control','is_sample'],'Smp':['is_sample']}

        if self.control_layout:
            # Parse LAYOUT ------------------------------------------------------
            if verbose > 0:
                logger.info(f"[TestPlate ApplyLayout] {self.plate_id} <- {self.control_layout} {self.n_wells}")
            _layLst = self.control_layout.split('_')
            _layDict = {}
            nLay = 0
            for _lo in CONTROL_ORDER:
                _l = _layLst[nLay]
                if _l != 'X':
                    _layDict[_lo] = {'R1':self.ROW_LABELS.index(_l[:1])+1,
                                     'C1':int(_l[1:3]),
                                     'R2':self.ROW_LABELS.index(_l[3:4])+1,
                                     'C2':int(_l[4:6])}
                else:
                    _layDict[_lo] = {}
                nLay += 1

            # ReSet LAYOUT ------------------------------------------------------
            for w in self.wells:
                for crt in CONTROL_LABELS:
                    self.set_well_field(w,crt,True)

            _n_layout = -1
            # Set per LAYOUT ------------------------------------------------------
            for _lo in CONTROL_ORDER:
                _rc = _layDict[_lo]
                if 'R1' in _rc :
                    for r in range(_rc['R1'],_rc['R2']+1):
                        for c in range(_rc['C1'],_rc['C2']+1):
                            for crt in CONTROL_ORDER[_lo]:
                                self.set_well_field((r,c),crt,True)

            self.n_layout= _n_layout 

    #--------------------------------------------------------------
    def calc_inhibition(self,verbose=0) -> int:
        if self.n_reads > 0:
            posReadOuts = self.get_readouts('is_poscontrol')
            pos_median = np.median(posReadOuts)
            pos_mad = np.median(np.absolute(posReadOuts - pos_median))
            pos_mean   = np.mean(posReadOuts)
            pos_std   = np.std(posReadOuts)

            negReadOuts = self.get_readouts('is_negcontrol')
            neg_median = np.median(negReadOuts)
            neg_mad = np.median(np.absolute(negReadOuts - neg_median))

            smpReadOuts = self.get_readouts('is_sample')
            smp_median = np.median(smpReadOuts)
            smp_mad    = np.median(np.absolute(smpReadOuts - smp_median))
            smp_mean   = np.mean(smpReadOuts)
            smp_std   = np.std(smpReadOuts)

            edgeReadOuts = self.get_readouts('is_edge')
            nonedgeReadOuts = self.get_readouts('is_nonedge')

            self.poscontrol_stats = [round(pos_median,4), round(pos_mad,4),round(pos_mean,4), round(pos_std,4)]
            self.negcontrol_stats = [round(neg_median,4), round(neg_mad,4),round(np.mean(negReadOuts),4), round(np.std(negReadOuts),4)]
            self.sample_stats     = [round(smp_median,4), round(smp_mad,4),round(smp_mean,4), round(smp_std,4)]
            self.edge_stats       = [np.median(edgeReadOuts), np.median(nonedgeReadOuts)]

            if verbose > 0:
                logger.info(f" {self.plate_id} Neg:{self.negcontrol_stats[0]} Pos:{self.poscontrol_stats[0]} Smp:{self.sample_stats[0]}")
                
            self.zfactor = round(1 - 3 * (pos_mad + neg_mad)/abs(pos_median - neg_median), 3)
            self.analysis_parameter = "Std pyAnalysis (dj)"

            _n_inhibition = 0
            for w in self.wells:
                self.wells[w].calc_inhibition(self.poscontrol_stats, self.negcontrol_stats,verbose=verbose)
                _n_inhibition += 1

            self.n_inhibition = _n_inhibition
            self.plate_qc = self.zfactor
            if self.test_issues:
                if 'Invalid' in self.test_issues:
                    self.plate_qc = -4
                if 'Dispensing' in self.test_issues:
                    self.plate_qc = -5
                if 'GrowthVariation' in self.test_issues:
                    self.plate_qc = -6
                if 'NoGrowth' in self.test_issues:
                    self.plate_qc = -6
                if 'Contamination' in self.test_issues:
                    self.plate_qc = -7

            if self.zfactor >= self.ZFACTOR_CUTOFF:
                self.plate_quality = setattr(self,'plate_quality',Dictionary.get(self.Choice_Dictionary['plate_quality'],'Valid')) 
            else:
                self.plate_quality = setattr(self,'plate_quality',Dictionary.get(self.Choice_Dictionary['plate_quality'],'Rejected'))
                self.test_issues = addto_StrList(self.test_issues,'FailedQC')

            if verbose > 0:
                _outstr  = f" {self.plate_quality} Zf: {self.zfactor:.3f} "
                _outstr += f"[POS: {self.poscontrol_stats[self.STATS_MEDIAN]:.3f} {self.poscontrol_stats[self.STATS_MAD]:.3f}] "
                _outstr += f"[NEG: {self.negcontrol_stats[self.STATS_MEDIAN]:.3f} {self.negcontrol_stats[self.STATS_MAD]:.3f}] "
                _outstr += f"[Edge: {self.edge_stats[0]:.3f} {self.edge_stats[1]:.3f}] "
                logger.info(f"[Calc Inhibition] {self.plate_id} - {_outstr} ")
        else:
            logger.warning(f"[Calc Inhibition] Plates has NO ReadOuts ")

    # -------------------------------------------------------
    def plot_heatmap(self,Property,outDir,propLegend=True):
    # -------------------------------------------------------
        
        if propLegend:
            bigTitle = f"{self.plate_id} ({self.result_type}) - {self.run_id} "
            subTitle = f"{Property} ({self.readout_type})"
            n_line = "\n"

            propTxt =  f"{n_line}Assay ID: {self.assay_id}"
            #propTxt += f"{n_line}Assay   : {self.PlateData['ASSAYTYPE_CODE']}"
            #propTxt += f"{n_line}Organism: {self.PlateData['ORGANISM']}"
            #propTxt += f"{n_line}Strain  : {self.PlateData['STRAIN']}"
            propTxt += f"{n_line}"
            propTxt += f"{n_line}Zfactor: {self.zfactor:.2f}"
            propTxt += f"{n_line}QC     : {self.plate_quality}"
            propTxt += f"{n_line}"
            propTxt += f"{n_line}PosCtrl: {self.poscontrol_stats[self.STATS_MEDIAN]:.2f}"
            propTxt += f"{n_line}NegCtrl: {self.poscontrol_stats[self.STATS_MEDIAN]:.2f}"
        else:
            bigTitle = f"{self.plate_id} (-) - {self.run_id} "
            subTitle = f"{Property} ({self.readout_type})"
            n_line = "\n"

            propTxt =  f"{n_line} - "

        #print(self.PlateData)
        gAxisY = f"Row"
        gAxisX = f"Column"

#=================================================================================================
class TestWell(Sample_Base):
    """

    """
#=================================================================================================

    Choice_Dictionary = {
        'conc_unit_lst':'Unit_Concentration',
        'conc_type_lst':'Concentration_Type',
        'solvent_conc_unit':'Unit_Concentration',
    }

    plate_id = models.ForeignKey(TestPlate, blank=False, verbose_name = "Plate ID", on_delete=models.DO_NOTHING,
        db_column="plate_id", related_name="%(class)s_plateid")
    well_id = models.CharField(max_length=5, blank=False, verbose_name = "Well ID")

    set_lst = ArrayField(models.CharField(max_length=5, blank=True),
                                 size=Sample_Base.MAX_CMPBATCHES, verbose_name = "Conc List", null=True, blank=True)

    volume = models.DecimalField(default=0, max_digits=12, decimal_places=4, verbose_name = "Volume [uL]")
    # volume_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Volume Unit", on_delete=models.DO_NOTHING,
    #      db_column="volume_unit", related_name="%(class)s_volume_unit")

    solvent = models.CharField(max_length=25, blank=True, verbose_name = "Solvent" )
    solvent_conc = models.DecimalField(default=0, max_digits=12, decimal_places=4, verbose_name = "SolvConc")
    solvent_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "SolvConc Unit", on_delete=models.DO_NOTHING,
         db_column="solvent_conc_unit", related_name="%(class)s_solvent_conc_unit")

    is_control = models.BooleanField(default=False, verbose_name = "is Control")
    is_poscontrol = models.BooleanField(default=False, verbose_name = "is PosCtrl")
    is_negcontrol = models.BooleanField(default=False, verbose_name = "is NegCtrl")
    is_sample = models.BooleanField(default=False, verbose_name = "is Sample")
    is_skip = models.BooleanField(default=False, verbose_name = "is Skip")
    is_valid = models.BooleanField(default=False, verbose_name = "is Valid")
    
    readouts = ArrayField(models.DecimalField(max_digits=12, decimal_places=5),
                          size=4, verbose_name = "Readouts", null=True, blank=True)
    readout_types = ArrayField(models.CharField(max_length=15, blank=True),
                          size=4, verbose_name = "Readout types", null=True, blank=True)
    
    inhibition = models.DecimalField(max_digits=9, decimal_places=3,default=-1)
    zscore = models.DecimalField(max_digits=9, decimal_places=3,default=-1)
    mscore = models.DecimalField(max_digits=9, decimal_places=3,default=-1)

    # I - Inactive, 
    # P - Partial (Inhib>=50 & MScore >= 2.5), 
    # A - Active  (Inhib>=80 & MScore >= 3.5)
    # S - Not Significant (Inhib>=80 & MScore <= 3.5)

    act_type = models.CharField(max_length=5, blank=True, verbose_name = "Act Type")

    # 0 - Inactive (I), 1 - Partial (P), 3 - Active (A) 
    act_score = models.SmallIntegerField(default=-1, blank=True, verbose_name = "Act Score")

    chk_migration = models.SmallIntegerField(default=-1, blank=False, verbose_name = "Check for migration")
    #-------------------------------------------------------------------------------
    class Meta:
        app_label = 'dplate'
        db_table = 'testwell'
        ordering=['plate_id','well_id']
        constraints = [
            models.UniqueConstraint(name='testwell_loc_cst', fields=['plate_id', 'well_id'], )
        ]        
        indexes = [
            GinIndex(name="testwell_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="testwell_wid_idx",fields=['well_id']),
            models.Index(name="testwell_is_idx",fields=['is_control', 'is_poscontrol', 'is_negcontrol', 'is_sample']),
            models.Index(name="testwell_skip_idx",fields=['is_skip', 'is_valid']),
            models.Index(name="testwell_atyp_idx",fields=['act_type']),
            models.Index(name="testwell_ascr_idx",fields=['act_score']),
            models.Index(name="testwell_inhib_idx",fields=['inhibition']),
            models.Index(name="testwell_mscr_idx",fields=['mscore']),
            models.Index(name="testwell_chkm_idx",fields=['chk_migration']),
        ]

    #-------------------------------------------------------------------------------
    def __str__(self):
        return f"{self.plate_id} {self.well_id}"

    #------------------------------------------------
    # Returns an TestWell instance if found by plate_id and well_id
    @classmethod
    def get(cls,PlateID,WellID,verbose=0):
        try:
            retInstance = cls.objects.get(plate_id=PlateID, well_id=WellID)
        except:
            if verbose:
                logger.warning(f"[Well Not Found] {PlateID} {WellID}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,PlateID,WellID):
        return cls.objects.filter(plate_id=PlateID, well_id=WellID).exists()

    # #------------------------------------------------
    def save(self, *args, **kwargs):
            if self.plate_id and self.well_id:
                super(TestWell, self).save(*args, **kwargs)
            else:
                logger.warning(f"[TestWell] SAVE has not PlateID and/or WellID") 

    #------------------------------------------------
    # 
    def get_welldict(self) -> dict:
        _Fields =[field.name for field in self._meta.fields if field.name not in self.AUDIT_FIELDS]
        return(model_to_dict(self, _Fields))

    #------------------------------------------------  
    def conv_list_to_string(self):
        super().conv_list_to_string()
        self.sets        = COMPOUND_SEP.join([str(x) for x in self.set_lst if x > 0])

    #------------------------------------------------  
    def conv_string_to_lst(self):
        super().conv_string_to_lst()
        self.set_lst = strList_to_List(self.sets,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def calc_inhibition(self,POS_Stats,NEG_Stats, verbose=0):
        _readout = float(self.readouts[0])
        _inhibition = 100 * (1 - (_readout - NEG_Stats[TestPlate.STATS_MEDIAN]) / (POS_Stats[TestPlate.STATS_MEDIAN] - NEG_Stats[TestPlate.STATS_MEDIAN]))
        # _zscore = (_readout - smp_mean) / smp_std
        # _mscore = 0.6745 * (_readout - smp_median) / smp_mad
        _zscore = (_readout - POS_Stats[TestPlate.STATS_MEAN]) / POS_Stats[TestPlate.STATS_STD]
        _mscore = 0.6745 * (_readout - POS_Stats[TestPlate.STATS_MEDIAN]) / POS_Stats[TestPlate.STATS_MAD]

        if verbose > 0:
            logger.info(f" {self.well_id} R:{_readout} N:{NEG_Stats} P:{POS_Stats} I:{_inhibition}")

        self.inhibition = round(_inhibition,2)
        self.zscore = round(_zscore,3)
        self.mscore = round(_mscore,3)

        _acttype = ActType_SC(_inhibition,_mscore)
        self.act_score = ActScoreSC_Cutoff[_acttype]['Score']
        self.act_type = ActScoreSC_Cutoff[_acttype]['Code']
