import os, re, math
import pandas as pd
import numpy as np
from decimal import Decimal

from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify
#from django.forms.models import model_to_dict

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from applib.data.str_lists import addto_StrList, strList_to_List
from dscreen.models import Screen_Run, Assay
from dorganism.models import Organism_Batch
from dcell.models import Cell_Batch
from dsample.models import Sample_Base, CmpBatchList_Base
from applib.bio.bio_data import ActScoreSC_Cutoff, ActType_SC
from adjcoadd.constants import *

import matplotlib.ticker as tic
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

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
    WELL_TYPE = Choices('Well','Tube','Vial')
    WELL_SIZE = Choices('Shallow','Deep','Storage')
    
    DICTIONARY_FIELDS = {
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
    
    PLATE_SIZES = {24:(4,6), 48:(6,8), 96:(8,12), 384:(16,24), 1536:(32,48), 400:(20,20), 2000:(25,80)}
    ROW_LABELS = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P',
                 'Q','R','S','T','U','V','W','X','Y','Z','AA','AB','AC','AD','AE','AF']
    MAP_POSITIONS = {'wellID':0,'pos2D':1,'pos1D':2}
    WELL_POS = MAP_POSITIONS['wellID']
    WELL_ID = MAP_POSITIONS['wellID']

    DICTIONARY_FIELDS = {
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
    def get(cls,PlateID, WellData=True, FillMissing=True, verbose=0):
        try:
            retInstance = cls.objects.get(plate_id=PlateID)
        except:
            if verbose:
                logger.warning(f"[Plate Not Found] {PlateID} ")
            retInstance = None

        if retInstance and WellData:
            retInstance.get_wells(FillMissing=FillMissing)

        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,PlateID,verbose=0):
        return cls.objects.filter(plate_id=PlateID).exists()

    #------------------------------------------------
    def set_platesize(self,PlateSize):
        # -- Set Plate Size
        if not np.isnan(PlateSize):
            if isinstance(PlateSize,float):
                PlateSize = int(PlateSize)

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
        _plate.plate_type = Dictionary.get(cls.DICTIONARY_FIELDS["plate_type"],PlateType)
        if WellData:
            _plate.init_wells(WellModel=None, PlateInstance=_plate)
        else:
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
            res = re.findall(r'([A-Za-z]+|\d+)',loc)
            loc = f"{res[0]}{int(res[1]):02d}"
            if check:
                if not loc in self.well_check['wellID']:
                    raise Exception(f"{self.plate_id} Invalid Well ID: {loc} for a {self.n_wells}w plate ({self.n_rows} x {self.n_cols})")
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
    def get_well(self,loc):
        if self.wells:
            w = self.map_well(loc)[self.WELL_POS]
            return(self.wells[w])

    #------------------------------------------------
    # Initialising Wells
    #------------------------------------------------
    def get_wells(self) -> int:
        # Create empty Wells
        self.init_wells()

        # Fill with Database Wells
        # to be implmented in specific models

    #--------------------------------------------------------------
    def load_wells(self, WellModel, FillMissing=True) -> int:
        # Create None Wells
        self.init_wells(WellModel=None, PlateInstance=None)

        # get Wells for that Plate 
        qryTW = WellModel.objects.filter(plate_id=self)
        lWells = qryTW.count()
        for w in qryTW:
            m = self.map_well(w.well_id)
            self.wells[m[0]] = w

        # Fill None Wells with empty {WellModel}
        if lWells < len(self.wells) and FillMissing:
            for w in self.wells:
                if self.wells[w] is None:
                    self.wells[w] = WellModel()
                    self.wells[w].well_id = w
                    self.wells[w].plate_id = self
            lWells = len(self.wells)
        
        return(lWells)

    #------------------------------------------------
    def save_wells(self) :
        if self.wells:
            for w in self.wells:
                if self.wells[w] is not None:
                    self.wells[w].save()
        
    #------------------------------------------------
    def delete_wells(self) :
        if self.wells:
            for w in self.wells:
                if self.wells[w]:
                    #print(f" {w} [{self.wells[w]}]")
                    self.wells[w].delete()
                    self.wells[w] = None

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
                    self.wells[m[self.WELL_POS]].well_id = m[self.WELL_ID]
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

    LIST_VIEW_FIELDS = {
        #"plate_id":{'Plate ID': {'plate_id':URL_LINKS['testplate_id']}},
        "plate_id":"Plate ID",
        "run_id":"Run ID",
        "assay_id":"Assay ID",
        "result_type":"Type",
        "control_layout":"Layout (N_P_R_S)",
        "test_date":"Test Date",
        "plate_quality":"Quality",
        "zfactor":"ZFactor",
        'n_reads' :"nR",
        'n_samples' : "nS",
        'n_layout' : "nL",
        'n_inhibitions' : "#SC",
        'n_doseresponses' : "#DR",
        'n_synergies'  : "#Syn",
        #"group_id":"Group",
        # "group_id.group_code":"Group",
    }

    WELL_CLASS = 'TestWell'
    
    RESULT_TYPES = Choices('MIC','CC50','HC50','SYN-MIC')
    ZFACTOR_CUTOFF = 0.2

    STATS_MEDIAN = 0
    STATS_MAD = 1
    STATS_MEAN = 2
    STATS_STD = 3

    DICTIONARY_FIELDS = {
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

    test_orgbatch_id = models.ForeignKey(Organism_Batch, null=True, blank=True, verbose_name = "OrgBatch", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatchid")
    test_cellbatch_id = models.ForeignKey(Cell_Batch, null=True, blank=True, verbose_name = "CellBatch", on_delete=models.DO_NOTHING,
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
    plate_comment = models.CharField(max_length=50, blank=True, verbose_name = "Plate Comment")
    
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
            _plate.plate_type = Dictionary.get(cls.DICTIONARY_FIELDS["plate_type"],'Test')

            if WellData:
                #print(f"[TestPlate.new] WithModel {TestWell}")
                _plate.init_wells(WellModel=TestWell,PlateInstance=_plate)
            else:
                _plate.init_wells(WellModel=None, PlateInstance=None)
            return(_plate)

    #------------------------------------------------
    def validate_model(self, WellData=True, **kwargs):
        retDict = []
        PlateDict = super(TestPlate, self).validate_model(**kwargs)
        for wd in PlateDict:
            retDict.append(wd)

        if hasattr(self,'wells') and WellData:
            for w in self.wells:
                if self.wells[w] is not None:
                    WellDict = super(TestWell,self.wells[w]).validate_model(**kwargs)
                    for wd in WellDict:
                        if 'plate_id' not in wd :
                            retDict.append(wd)
        
        return(retDict)

    #------------------------------------------------
    def set_defaults_model(self, WellData=True, verbose = 0):
        retDict = []
        super(TestPlate, self).set_defaults_model()

        # if self.wells have been set at all
        if hasattr(self,'wells'):
            # if self.wells is not None 
            if self.wells and WellData:
                for w in self.wells:
                    if self.wells[w] is not None:
                        super(TestWell,self.wells[w]).set_defaults_model()
        
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
    def get_wells(self, FillMissing=True) -> int:
        
        # Create None Wells
        self.init_wells(WellModel=None, PlateInstance=None)

        # get Wells for that Plate 
        qryTW = TestWell.objects.filter(plate_id=self)
        lWells = qryTW.count()
        for w in qryTW:
            m = self.map_well(w.well_id)
            self.wells[m[0]] = w

        # Fill None Wells with empty TestWell
        if lWells < len(self.wells) and FillMissing:
            for w in self.wells:
                if self.wells[w] is None:
                    self.wells[w] = TestWell()
                    self.wells[w].well_id = w
                    self.wells[w].plate_id = self
            lWells = len(self.wells)
        
        return(lWells)

    #--------------------------------------------------------------
    def make_wells_df(self, RowCol=False, ListToString=False, ReadoutField=True) -> pd.DataFrame:
        _dicts = []
        if self.wells:
            for w in self.wells:
                if self.wells[w] is not None:
                    _well_dict = self.wells[w].get_well_dict(ListToString=ListToString, ReadoutField=ReadoutField)
                    if RowCol:
                        _r,_c = self.well_rowcol(w)
                        _well_dict['row'] = self.ROW_LABELS[_r-1]
                        _well_dict['col'] = _c
        
                    _dicts.append(_well_dict)
        self.wells_df = pd.DataFrame(_dicts)
        return(len(self.wells_df))

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
    def update_n(self,nType='n_samples'):
        if nType == 'n_samples':
            n_sample = 0
            if self.n_wells >0:
                for w in self.wells:
                    if self.wells[w].n_cmpbatches > 0:
                        n_sample += 1
            self.n_samples = n_sample

    #--------------------------------------------------------------
    def clear_cmpbatch_data(self):  
        #
        # resets cmpbatch data (incl conc, conc_unit, conc_type)
        if hasattr(self,'wells'):
            for w in self.wells:
                self.wells[w].clear_cmpbatch_data()

    #--------------------------------------------------------------
    def conv_list_to_string(self):
        if hasattr(self,'wells'):
            for w in self.wells:
                self.wells[w].conv_list_to_string()
        
    #--------------------------------------------------------------
    def process_wells(self,Function):
        if hasattr(self,'wells'):
            for w in self.wells:
                self.wells[w].Function()
 
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

        _n_layout = -1
        if self.control_layout:
            # Parse LAYOUT ------------------------------------------------------
            if verbose > 0:
                logger.info(f"[TestPlate ApplyLayout] {self.plate_id} <- {self.control_layout} {self.n_wells}")
            _layLst = self.control_layout.split('_')
            _layDict = {}
            nLay = 0
            for _lo in CONTROL_ORDER:
                _l = _layLst[nLay]
#                if _l != 'X' or _l != 'MIC':
                if _l not in ['X','MIC']:
                    _layDict[_lo] = {'R1':self.ROW_LABELS.index(_l[:1])+1,
                                     'C1':int(_l[1:3]),
                                     'R2':self.ROW_LABELS.index(_l[3:4])+1,
                                     'C2':int(_l[4:6])}
                else:
                    _layDict[_lo] = {}
                nLay += 1

            # ReSet LAYOUT ------------------------------------------------------
            _n_layout = -1
            for w in self.wells:
                for crt in CONTROL_LABELS:
                    self.set_well_field(w,crt,False)
            # Set per LAYOUT ------------------------------------------------------
            _n_layout = 0
            for _lo in CONTROL_ORDER:
                _rc = _layDict[_lo]
                if 'R1' in _rc :
                    _n_layout += 1
                    for r in range(_rc['R1'],_rc['R2']+1):
                        for c in range(_rc['C1'],_rc['C2']+1):
                            for crt in CONTROL_ORDER[_lo]:
                                self.set_well_field((r,c),crt,True)

        self.n_layout= _n_layout
        return(_n_layout) 

    #--------------------------------------------------------------
    def calc_inhibition(self,verbose=0) -> int:
        _nInhibs = 0
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

            
            for w in self.wells:
                self.wells[w].calc_inhibition(self.poscontrol_stats, self.negcontrol_stats,verbose=verbose)
                _nInhibs += 1

            self.n_inhibitions = _nInhibs

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
                setattr(self,'plate_quality',Dictionary.get(self.DICTIONARY_FIELDS['plate_quality'],'Valid')) 
            else:
                setattr(self,'plate_quality',Dictionary.get(self.DICTIONARY_FIELDS['plate_quality'],'Rejected'))

                self.test_issues = addto_StrList(self.test_issues,'FailedQC')

            if verbose > 0:
                _outstr  = f" {self.plate_quality} Zf: {self.zfactor:.3f} "
                _outstr += f"[POS: {self.poscontrol_stats[self.STATS_MEDIAN]:.3f} {self.poscontrol_stats[self.STATS_MAD]:.3f}] "
                _outstr += f"[NEG: {self.negcontrol_stats[self.STATS_MEDIAN]:.3f} {self.negcontrol_stats[self.STATS_MAD]:.3f}] "
                _outstr += f"[Edge: {self.edge_stats[0]:.3f} {self.edge_stats[1]:.3f}] "
                logger.info(f"[Calc Inhibition] {self.plate_id} - {_outstr} ")
        else:
            logger.warning(f"[Calc Inhibition] Plates has NO ReadOuts ")

        return(_nInhibs)

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
            propTxt += f"{n_line}NegCtrl: {self.negcontrol_stats[self.STATS_MEDIAN]:.2f}"
        else:
            bigTitle = f"{self.plate_id} (-) - {self.run_id} "
            subTitle = f"{Property} ({self.readout_type})"
            n_line = "\n"

            propTxt =  f"{n_line} - "

        #print(self.PlateData)
        gAxisY = f"Row"
        gAxisX = f"Column"

        if not hasattr(self, 'wells_df'):
            self.make_wells_df(RowCol=True)

        self.wells_df = self.wells_df.astype({Property: 'float'})
        prop_map = self.wells_df.pivot_table(index="row", columns="col", values=Property)

        fig, ax = plt.subplots(figsize=(12,6))
        fig.text(0.05,0.91,bigTitle, fontsize=19, ha = 'left')
        fig.text(0.97,0.91,subTitle, fontsize=12, color = 'darkgrey', ha = 'right')
        fig.text(0.84,0.80,propTxt,fontsize=9, color = 'black', ha = 'left', va='top',wrap=True,backgroundcolor='lightgray')
        
        if Property == 'inhibition':
            _fmt = ".1f"
            _vmin = 0
            _vmax = 100
            _col = plt.cm.get_cmap('RdYlGn_r')
        else:
            _fmt = ".3f"
            #_vmin,_vmax = well_df[Property].quantile([.01, .99])
            _vmax = self.wells_df[Property].max()
            _vmin = self.wells_df[Property].min()
            _col = sns.light_palette("darkred", as_cmap=True)
        
        ax = sns.heatmap(prop_map, 
                         linewidth=.5, 
                         annot=True, fmt=_fmt, annot_kws={'size': 5},
                         xticklabels=True, yticklabels=True,
                         vmin=_vmin, vmax=_vmax, 
                         cmap=_col)
        
        plt.yticks(rotation=0) 
        plt.xlabel(gAxisX, fontsize= 12)
        plt.ylabel(gAxisY, fontsize= 12)

        if outDir:
            xOutDir = os.path.join(outDir,str(self.run_id))
            if not os.path.exists(xOutDir):
                os.makedirs(xOutDir)

            jpgFile = f"{self.plate_id}_{Property}.jpg"
            fig.savefig(os.path.join(xOutDir,jpgFile))
            plt.close(fig)
        else:
            fig.show()


#=================================================================================================
class TestWell(Sample_Base):
    """

    """
#=================================================================================================

    DICTIONARY_FIELDS = {
        'conc_unit_lst':'Unit_Concentration',
        'conc_type_lst':'Concentration_Type',
        'solvent_conc_unit':'Unit_Concentration',
    }

    STRING_FIELDS = Sample_Base.STRING_FIELDS + ['sets','cmpbatch_sets']

    ARRAY_FIELDS = {'cmpbatch_lst':['compound_id','compound2_id','compound3_id','compound4_id'],
                    'conc_lst':['conc','conc2','conc3','conc4',],
                    'conc_unit_lst':['conc_unit','conc2_unit','conc3_unit','conc4_unit'], 
                    'conc_type_lst':['conc_type','conc2_type','conc3_type','conc4_type'], 
                    'set_lst':['set_id','set2_id','set3_id','set4_id'], 
                    }
    
    COPY_FIELDS = ['solvent', 'solvent_conc', 'amount','volume',]

    # Fields ---------------------------------------------------------------------------------------------------

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
    
    # CmpBatch+Set list - for Doseresponse
    cmpbatch_sets = ''
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

    def str_cmpbatch_data(self):
        return f"{self.plate_id} {self.well_id} {self.cmpbatch_lst} {self.conc_lst} {self.conc_unit_lst}"

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
    def conv_list_to_string(self):
        super().conv_list_to_string()

        self.sets = ''
        if self.set_lst:
            self.sets = COMPOUND_SEP.join([str(x) for x in self.set_lst])
        self.cmpbatch_sets = ''
        if self.set_lst:
            self.cmpbatch_sets = COMPOUND_SEP.join([f"{str(c)}_{str(s)}" for c,s in zip(self.cmpbatch_lst, self.set_lst) ])
        else:
            self.cmpbatch_sets = self.cmpbatches
            
    #------------------------------------------------  
    def conv_string_to_lst(self):
        super().conv_string_to_lst()
        self.set_lst = strList_to_List(self.sets,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def fields_to_dict(self,AuditFields=False, ModelFields=[], ClassFields=[], ReadoutField=True):
        _dict = super().fields_to_dict(AuditFields=AuditFields, ModelFields=ModelFields, ClassFields=ClassFields)
        if ReadoutField:
            if self.readouts:
                for i in range(len(self.readouts)):
                    _dict[f'readout_{i+1}'] = self.readouts[i]
        return(_dict)

    #------------------------------------------------
    def get_well_dict(self, ListToString=False, ReadoutField=True) -> dict:
        _ClassFields = []
        if ListToString:
            self.conv_list_to_string()
            _ClassFields += self.STRING_FIELDS
        return(self.fields_to_dict(ClassFields=_ClassFields,ReadoutField=ReadoutField))
    
    #------------------------------------------------
    def clear_cmpbatch_data(self):
        super().clear_cmpbatch_data()
        self.set_lst = []
    
    #------------------------------------------------
    def clear_inhibition_data(self):
        self.inhibition = -1
        self.zscore = -1
        self.mscore = -1
        self.act_score = -1
        self.act_type = ''
            

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

#=================================================================================================#=================================================================================================
class MasterPlate(Plate):
    """

    """
#=================================================================================================

    DICTIONARY_FIELDS = {
        'plate_quality':'Data_Quality',
        'plate_type':'Plate_Type',
        'master_conc_unit':'Unit_Concentration',
        'master_volume_unit':'Unit_Volume',
    }

    plating = models.CharField(max_length=10, blank=True, verbose_name = "Plating by")
    cpoz_id = models.CharField(max_length=15, blank=True, verbose_name = "CpOz ID")
    well_type= models.CharField(max_length=20, blank=True, null=True, choices=Labware.WELL_TYPE, verbose_name = "Type")  

    run_id = models.ForeignKey(Screen_Run, null=True, blank=True, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_runid")    
    prep_date = models.DateField(null=True, blank=True, verbose_name = "Prep Date")

    plate_quality = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Quality", on_delete=models.DO_NOTHING,
        db_column="plate_quality", related_name="%(class)s_platequality")
    plate_comment = models.CharField(max_length=50, blank=True, verbose_name = "Plate Comment")

#    control_layout = models.CharField(max_length=35, blank=True, verbose_name = "Control Layout")
    dilution_layout = models.CharField(max_length=35, blank=True, verbose_name = "Dilution Layout")

    n_samples = models.SmallIntegerField(default=-1, verbose_name = "#Samples")
    n_layout = models.SmallIntegerField(default=-1, verbose_name = "Layout")
    process_status = models.SmallIntegerField(default=-1, verbose_name = "Process Status")

    # master_processing = models.CharField(max_length=25, blank=True, verbose_name = "Processing")
    # master_issues = models.CharField(max_length=150, blank=True, verbose_name = "Issue")

    master_solvent = models.CharField(max_length=25, blank=True, verbose_name = "Solvent" )
    master_conc = models.DecimalField(default=-1, max_digits=12, decimal_places=4, verbose_name = "Conc")
    master_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Conc Unit", on_delete=models.DO_NOTHING,
         db_column="master_conc_unit", related_name="%(class)s_master_conc_unit")
    master_volume = models.DecimalField(default=-1, max_digits=10, decimal_places=2,  verbose_name = "Volume (uL)")
    master_volume_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Volume Unit", on_delete=models.DO_NOTHING,
         db_column="master_volume_unit", related_name="%(class)s_master_volume_unit")

    class Meta:
        app_label = 'dplate'
        db_table = 'masterplate'
        ordering=['plate_id']
        indexes = [
            models.Index(name="masterplate_labw_idx",fields=['labware_id']),
            models.Index(name="masterplate_wtyp_idx",fields=['well_type']),
            models.Index(name="masterplate_run_idx",fields=['run_id']),
            models.Index(name="masterplate_proc_idx",fields=['process_status']),
            models.Index(name="masterplate_pq_idx",fields=['plate_quality']),
        ]

    #--------------------------------------------------------------
    def __repr__(self):
        _str  = f" [Masterplate] {self.plate_id} Size:{self.n_wells} {self.plate_type}"
        if hasattr(self,'wells'):
            _wellid = list(self.wells.keys())
            _str += f" Wells:{len(self.wells)} [{_wellid[0]}..{_wellid[-1]}]"           
        return(_str)


    #------------------------------------------------
    @classmethod
    def new(cls,PlateID,PlateSize,PlateType,WellData=True,NoCheck=False):
        #print(f"[MasterPlate.new] {PlateID} {PlateSize} ")
        if not cls.exists(PlateID.upper()) or NoCheck:
            _plate = cls()
            _plate.plate_id = PlateID.upper()
            _plate.set_platesize(PlateSize)
            _plate.plate_type = Dictionary.get(cls.DICTIONARY_FIELDS["plate_type"],PlateType)

            if WellData:
                #print(f"[MasterPlate.new] WithModel {TestWell}")
                _plate.init_wells(WellModel=MasterWell,PlateInstance=_plate)
            else:
                _plate.init_wells(WellModel=None, PlateInstance=None)
            return(_plate)
        else:
            logger.warning(f"[Masterplate] New {PlateID.upper()} alreday exists ")
            return(None)

    #--------------------------------------------------------------
    def get_wells(self, FillMissing=True) -> int:
        # Create None Wells
        self.init_wells(WellModel=None, PlateInstance=None)

        # get Wells for that Plate 
        qryTW = MasterWell.objects.filter(plate_id=self)
        lWells = qryTW.count()
        for w in qryTW:
            m = self.map_well(w.well_id)
            self.wells[m[0]] = w

        # Fill None Wells with empty TestWell
        if lWells < len(self.wells) and FillMissing:
            for w in self.wells:
                if self.wells[w] is None:
                    self.wells[w] = MasterWell()
                    self.wells[w].well_id = w
                    self.wells[w].plate_id = self
            lWells = len(self.wells)
        
        return(lWells)

    #------------------------------------------------
    def validate_model(self, WellData=True, **kwargs):
        retDict = []
        PlateDict = super(MasterPlate, self).validate_model(**kwargs)
        for wd in PlateDict:
            retDict.append(wd)

        if hasattr(self,'wells') and WellData:
            for w in self.wells:
                if self.wells[w] is not None:
                    WellDict = super(MasterWell,self.wells[w]).validate_model(**kwargs)
                    for wd in WellDict:
                        if 'plate_id' not in wd :
                            retDict.append(wd)
        
        return(retDict)

    #------------------------------------------------
    def set_defaults_model(self, WellData=True, verbose = 0):
        retDict = []
        super(MasterPlate, self).set_defaults_model()

        if hasattr(self,'wells') and WellData:
            for w in self.wells:
                if self.wells[w] is not None:
                    super(MasterWell,self.wells[w]).set_defaults_model(ignore_fields=['barcode'])
        
    #------------------------------------------------
    def save(self, *args, **kwargs):
        if self.plate_id:
            verbose = kwargs.get('verbose',0)
            kwargs.pop("verbose",None)

            super(MasterPlate, self).save(*args, **kwargs)
            if verbose > 0:
                logger.info(f"[MasterPlate.save] {self.plate_id}")
            if hasattr(self,'wells'):
                for w in self.wells:
                    if self.wells[w] is not None:
                        if verbose > 1:
                            logger.info(f"[MasterWell.save] {self.wells[w]}")        
                        super(MasterWell,self.wells[w]).save(*args, **kwargs)
        else:
            logger.warning(f"[MasterPlate] SAVE has no PlateID ")

    #------------------------------------------------
    def add_dilutions(self, **kwargs):
        DILUTION_DICT = {
            'Col8':     ( 8, 2, True, False),
            'Col16':    (16, 2, True, False),
            'Fix_Col16':(16, 0, True, False),
            'Fix_Col8': ( 8, 0, True, False),
            'Row8':     ( 8, 2, False, True),
            'Row10':    (10, 2, False, True),
        }

        verbose = kwargs.get('verbose',0)
        valLog = kwargs.get('valLog',None)
       
        if self.wells:
            for w in self.wells:
                if self.wells[w].dilution_lst:
                    if verbose > 0:
                        logger.info(f"[MasterPlate] Dilution [{self.plate_id} {w}] {self.wells[w].dilution_lst} {self.wells[w].test_conc_lst}")

                    for i_dil in range(len(self.wells[w].dilution_lst)):

                        d =self.wells[w].dilution_lst[i_dil]

                        if d in DILUTION_DICT:
                            nConc,dConc,dRow,dCol  = DILUTION_DICT[d]
                            w_R,w_C = self.well_rowcol(w)
                            w_Conc = self.wells[w].test_conc_lst[i_dil]

                            #print(f"{self.wells[w].plate_id} {self.wells[w].well_id} {self.wells[w].test_conc_lst} [{i_dil}] {d} -> {w_Conc}  ")

                            for n in range(nConc-1):
                                if dConc > 0:
                                    w_Conc = w_Conc / dConc
                                if dRow:
                                    w_R += 1
                                elif dCol:
                                    w_C += 1
                                w_d = self.well_id((w_R,w_C))
                                
                                # Set Dilution Well if empty
                                if self.wells[w_d].n_cmpbatches == 0:
                                    #print(f"{self.plate_id} - {w_d} {i_dil} Reset")
                                    self.wells[w_d].cmpbatch_lst = self.wells[w].cmpbatch_lst
                                    self.wells[w_d].n_cmpbatches = self.wells[w].n_cmpbatches
                                    _conc = list(getattr(self.wells[w],'test_conc_lst'))
                                    setattr(self.wells[w_d],'test_conc_lst',_conc)
                                    #self.wells[w_d].test_conc_lst = self.wells[w].test_conc_lst
                                    self.wells[w_d].test_conc_unit_lst = self.wells[w].test_conc_unit_lst
                                    self.wells[w_d].set_lst = self.wells[w].set_lst
                                    #self.wells[w_d].dilution_lst = []
                                    
                                # Set test_conc of i-th cmpbatch to wconc
                                # _wconc_lst = getattr(self.wells[w_d],'test_conc_lst')
                                # _wconc_lst[i_dil] = w_Conc
                                # print(_wconc_lst)
                                # setattr(,'test_conc_lst',_wconc_lst)
                                self.wells[w_d].test_conc_lst[i_dil] = w_Conc

                                #print(f"{self.wells[w_d].plate_id} {self.wells[w_d].well_id} : [{w_d}] {self.wells[w_d].test_conc_lst} <--  [{w}] {self.wells[w].test_conc_lst} ")

                                # for x in ['A19','B19','C19','D19']:
                                #     print(f" xxx  {self.plate_id}  {self.wells[x].well_id} {self.wells[x].test_conc_lst}")

                                #print(f" {i} {self.plate_id} {dw_id} {wConc} {self.wells[dw_id].test_conc_lst}")
                        else:
                            logger.warning(f"[MasterPlate] Unknown dilution {d} [{self.plate_id} {w}]")
                            if valLog:
                                valLog.add_error("Unknown dilution",f"Diluation [{d}]",f"{self.plate_id}:{w}",f"Correct Dilution {list(DILUTION_DICT.keys())}")

                    


#=================================================================================================
class MasterWell(Sample_Base):
    """

    """
#=================================================================================================

    #STRING_FIELDS =['concs','conc_units','conc_types','sets']
    STRING_FIELDS = Sample_Base.STRING_FIELDS + ['sets','test_concs','test_conc_units','test_conc_types']

    DICTIONARY_FIELDS = {
        'conc_unit_lst':'Unit_Concentration',
        'conc_type_lst':'Concentration_Type',
        'solvent_conc_unit':'Unit_Concentration',
        'amount_unit':'Unit_Amount',
        'volume_unit':'Unit_Volume',
        'test_conc_unit_lst': 'Unit_Concentration',
        'test_conc_type_lst': 'Concentration_Type',
    }

    ARRAY_FIELDS = {'cmpbatch_lst':['compound_id','compound2_id','compound3_id','compound4_id'],
                    'conc_lst':['conc','conc2','conc3','conc4',],
                    #'conc_unit_lst':['conc_unit','conc2_unit','conc3_unit','conc4_unit'], 
                    #'conc_type_lst':['conc_type','conc2_type','conc3_type','conc4_type'], 
                    'set_lst':['set_id','set2_id','set3_id','set4_id'], 
                    'test_conc_lst':['test_conc','test_conc2','test_conc3','test_conc4',],
                    #'test_conc_unit_lst':['test_conc_unit','test_conc2_unit','test_conc3_unit','test_conc4_unit'], 
                    #'test_conc_type_lst':['test_conc_type','test_conc2_type','test_conc3_type','test_conc4_type'], 
                    'dilution_lst':['dilution','dilution2','dilution3','dilution4'], 
                    }
    
    ARRAYDICTIONARY_FIELDS = {
                    'conc_unit_lst':['conc_unit','conc2_unit','conc3_unit','conc4_unit'], 
                    'conc_type_lst':['conc_type','conc2_type','conc3_type','conc4_type'], 
                    'test_conc_unit_lst':['test_conc_unit','test_conc2_unit','test_conc3_unit','test_conc4_unit'],     
                    'test_conc_type_lst':['test_conc_type','test_conc2_type','test_conc3_type','test_conc4_type'],      
    }
    
    COPY_FIELDS = ['solvent', 'solvent_conc', 'amount','volume',]


    plate_id = models.ForeignKey(MasterPlate, blank=True, null=True, verbose_name = "Plate ID", on_delete=models.DO_NOTHING,
        db_column="plate_id", related_name="%(class)s_plateid")
    well_id = models.CharField(max_length=5, blank=True, null=True, verbose_name = "Well ID")
    barcode = models.CharField(max_length=15, blank=True, null=True, unique=True, verbose_name = "Barcode")

    sets = ""
    set_lst = ArrayField(models.CharField(max_length=5, blank=True),
                                 size=Sample_Base.MAX_CMPBATCHES, verbose_name = "Conc List", null=True, blank=True)

    solvent = models.CharField(max_length=50, blank=True, verbose_name = "Solvent" )
    solvent_conc = models.DecimalField(default=0, max_digits=12, decimal_places=4, verbose_name = "SolvConc")
    solvent_conc_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "SolvConc Unit", on_delete=models.DO_NOTHING,
         db_column="solvent_conc_unit", related_name="%(class)s_solvent_conc_unit")

    amount = models.DecimalField(default=-1, max_digits=12, decimal_places=4, verbose_name = "Amount")
    amount_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Amount Unit", on_delete=models.DO_NOTHING,
         db_column="amount_unit", related_name="%(class)s_amount_unit")
    volume = models.DecimalField(default=-1, max_digits=10, decimal_places=2,  verbose_name = "Volume")
    volume_unit = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Volume Unit", on_delete=models.DO_NOTHING,
         db_column="volume_unit", related_name="%(class)s_volume_unit")

    dilution_lst = ArrayField(models.CharField(max_length=15, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "Dilution List", null=True, blank=True)

    test_concs = ""
    test_conc_lst = ArrayField(models.DecimalField(max_digits=9, decimal_places=4, default=0), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "TestConc List", null=True, blank=True)
    test_conc_units = ""
    test_conc_unit_lst = ArrayField(models.CharField(max_length=10, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "TestConcUnit List", null=True, blank=True)
    test_conc_types = ""
    test_conc_type_lst = ArrayField(models.CharField(max_length=5, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "TestConcType List", null=True, blank=True)

    prev_plate_id = models.CharField(max_length=25, blank=True, verbose_name = "Prev Plate ID")
    prev_well_id = models.CharField(max_length=5, blank=True, verbose_name = "Prev Well ID")

    chk_migration = models.SmallIntegerField(default=-1, blank=False, verbose_name = "Check for migration")

# CPOZ2_SN	VARCHAR2(25 BYTE)
# CPOZ_SN	VARCHAR2(25 BYTE)
# TEST_SOLVENT_CONC	NUMBER

    class Meta:
        app_label = 'dplate'
        db_table = 'masterwell'
        ordering=['plate_id','well_id']
        constraints = [
            models.UniqueConstraint(name='masterwell_loc_cst', fields=['plate_id', 'well_id'], )
        ]        
        indexes = [
            GinIndex(name="masterwell_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="masterwell_wid_idx",fields=['well_id']),
            models.Index(name="masterwell_bc_idx",fields=['barcode']),
            #models.Index(name="testwell_is_idx",fields=['is_control', 'is_poscontrol', 'is_negcontrol', 'is_sample']),
            models.Index(name="masterwell_chkm_idx",fields=['chk_migration']),
        ]
    #-------------------------------------------------------------------------------
    def __str__(self):
        return f"{self.plate_id} {self.well_id} {self.barcode}"
    
    def str_cmpbatch_data(self):
        return f"{self.plate_id} {self.well_id} {self.cmpbatch_lst} {self.test_conc_lst} {self.test_conc_unit_lst}"

    #------------------------------------------------
    # Returns an MasterWell instance if found by plate_id and well_id
    @classmethod
    def get(cls,PlateID,WellID,Barcode=None,verbose=0):
        try:
            if Barcode:
                retInstance = cls.objects.get(barcode=Barcode)
            else:
                retInstance = cls.objects.get(plate_id=PlateID, well_id=WellID)
        except:
            if verbose:
                logger.warning(f"[Well Not Found] {PlateID} {WellID} {Barcode}")
            retInstance = None
        return(retInstance)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,PlateID,WellID,Barcode=None):
        if Barcode:
            return cls.objects.filter(barcode=Barcode).exists()
        else:    
            return cls.objects.filter(plate_id=PlateID, well_id=WellID).exists()

    # #------------------------------------------------
    def save(self, *args, **kwargs):
            if (self.plate_id and self.well_id) or self.barcode:
                
                # Sets empty barcodes ('') to None 
                #    required for unique barcode constrain 
                self.set_none_field('barcode')
                
                super(MasterWell, self).save(*args, **kwargs)
            else:
                logger.warning(f"[MasterWell] SAVE has not PlateID and/or WellID") 


    #------------------------------------------------  
    def conv_list_to_string(self):
        super().conv_list_to_string()

        self.sets = ''
        if self.set_lst:
            self.sets = COMPOUND_SEP.join([str(x) for x in self.set_lst != ""])

        self.test_concs = ''
        if self.test_conc_lst:
            self.test_concs  = COMPOUND_SEP.join([str(x) for x in self.test_conc_lst if x > 0])
        self.test_conc_units = ''
        if self.test_conc_unit_lst:
            self.test_conc_units = COMPOUND_SEP.join([str(x) for x in self.test_conc_unit_lst != ""])
        self.test_conc_types = ''
        if self.test_conc_type_lst:
            self.test_conc_types = COMPOUND_SEP.join([str(x) for x in self.test_conc_type_lst != ""])

    #------------------------------------------------  
    def conv_string_to_lst(self):
        super().conv_string_to_lst()
        self.set_lst = strList_to_List(self.sets,sep=COMPOUND_SEP,size=4,fill="")
        self.test_conc_lst = strList_to_List(self.test_concs,sep=COMPOUND_SEP,size=4,fill=0)
        self.test_conc_units = strList_to_List(self.test_conc_units,sep=COMPOUND_SEP,size=4,fill="")
        self.test_conc_types = strList_to_List(self.test_conc_types,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------
    def clear_cmpbatch_data(self):
        super().clear_cmpbatch_data()
        self.set_lst = []
        self.test_concs = ""
        self.test_conc_lst = []
        self.test_conc_units = ""
        self.test_conc_unit_lst = []
        self.test_conc_types = ""
        self.test_conc_type_lst = []
