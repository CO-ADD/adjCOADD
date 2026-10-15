import sys, os
import datetime
import numpy as np
import pandas as pd
import scipy.optimize as opt

from adjcoadd.constants import COMPOUND_SEP
from apputil.models import Dictionary
from applib.bio.bio_data import ActType_DR, pScore, format_DR, dr_max_quality, ActScoreDR_Cutoff
from dsample.models import Compound_Batch
from dplate.models import TestPlate
from dscreen.models import AssayData_MIC,AssayData_CC50,AssayData_HC50
from adjcoadd.constants import DR_CLASSES
import logging
logger = logging.getLogger(__name__)

#====================================================================
class Synergy():
#====================================================================
    
    #--------------------------------------------------------------
    def __init__(self, inhibition_cutoff=80, 
                 min_dilutions=6,  inhib_limit = 500, 
                 inhib_correction = True, inhib_correct_limit = 20,
                 ic50_dmax_cutoff = 40,
                 breakpoint=None, 
                 ):
        self.dr_quality = 'Empty'
        self.dr_type = None
        self.dr_dmax = ''

        self.inhibition_cutoff = inhibition_cutoff
        self.min_dilutions = min_dilutions
        self.breakpoint = breakpoint
        self.inhib_limit = inhib_limit
        self.inhib_correct = inhib_correction 
        self.inhib_correct_limit = inhib_correct_limit
        self.ic50_dmax_cutoff = ic50_dmax_cutoff
