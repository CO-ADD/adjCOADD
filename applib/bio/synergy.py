import sys, os
import datetime
import numpy as np
import pandas as pd
import scipy.optimize as opt

from adjcoadd.constants import COMPOUND_SEP
from apputil.models import Dictionary
from applib.bio.bio_data import ActType_DR, pScore, format_DR, dr_max_quality, ActScoreDR_Cutoff
from applib.bio.doseresponse import DoseResponse
from dsample.models import Compound_Batch
from dplate.models import TestPlate, TestWell
from dscreen.models import AssayData_MIC,AssayData_CC50,AssayData_HC50
from adjcoadd.constants import SYN_CLASSES

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

    #--------------------------------------------------------------
    # Initialise Synergy with
    #   CmpBatches - string of cmpbatches
    #   dfDR - dataframe of dosereponse data
    #   djTestPlate - TestPlate instance

    @staticmethod
    def apply_conc(s):
        if s['conc_lst']:
            s['conc'] = s['conc_lst'][0]
        else:
            logger.warning(" [Synergy] No concentrations ")
        if s['conc_unit_lst']:
            s['conc_unit'] = s['conc_unit_lst'][0]
        else:
            logger.warning(" [Synergy] No concentration unit ")
        return(s)

    def init_data(self,CmpBatches, dfSyn, djTestPlate):
        dfValid = True
        self.required_columns = ['well_id','conc_lst', 'inhibition','conc_unit_lst']
        # Prepare Data DF
        if set(self.required_columns).issubset(dfSyn.columns) :
            dfValid = True


# Process TestPlate for Synergies ------------------------------------------------
# 
def process_testplate_synergy(djTP,upload=False,overwrite=False,verbose=0):
    DR_Lst = []
    #djTP = TestPlate.get(PlateID,WellData=True)
    if djTP:
        if str(djTP.result_type) in SYN_CLASSES:
            if djTP.n_samples > 0 and djTP.n_inhibitions > 0 :
                djTP.n_synergy = 0
                if verbose > 0:
                    logger.info(f" [{djTP.plate_id}] {djTP.result_type} {djTP.assay_id} ")
                djTP.conv_list_to_string()
                djTP.make_wells_df(ListToString=True)

                _synDrugs = []
                _synPots = []
                if djTP.synergy_cmpbatches[0]:
                    _synDrugs = djTP.synergy_cmpbatches[0].split(';')
                if djTP.synergy_cmpbatches[1]:
                    _synPots = djTP.synergy_cmpbatches[1].split(';')
                logger.debug(f" [{djTP.plate_id}] --1-- Drugs: {_synDrugs} - Potentiators: {_synPots}")
                
                print(djTP.wells_df)
                grpData = djTP.wells_df.groupby('cmpbatch_sets')
                
                for CmpBatchSet,SynData in grpData:
                    if CmpBatchSet:
                        # Single MIC for Drug and Potentiator
                        #--------------------------------------
                        _cmpdsets = CmpBatchSet.split(COMPOUND_SEP)
                        if len(_cmpdsets) == 1:
                            _cmpbatch = TestWell.split_cmpbatch_sets(_cmpdsets[0])[0]
                            djDR = DoseResponse().init_data(_cmpbatch,SynData[['well_id','conc_lst','inhibition','mscore','conc_unit_lst']],djTP)
                            if djDR:
                                djDR.calc_doseresponse()
                                print(repr(djDR))
                                djTP.n_doseresponses += 1
                                
                                djDR.doseresponse_to_assaydata() 
                                #print(f" {djDR.assaydata.data_quality} ({djDR.assaydata.data_comment}) ")
                                if upload:
                                    djDR.save_assaydata(overwrite=overwrite)
                                DR_Lst.append(djDR)

                for CmpBatchSet,SynData in grpData:
                    if CmpBatchSet:
                        # Synergy MIC for Drug|Potentiator
                        #--------------------------------------
                        _cmpdsets = CmpBatchSet.split(COMPOUND_SEP)
                        if len(_cmpdsets) > 1:
                            print(_cmpdsets, _synDrugs, _synPots)
                            
                        #     _drug, _pot = _cmpdsets
                        #     #
                        # else:
                        #     _cmpd = TestWell.split_cmpbatch_sets(_cmpdsets[0])[0]
                        #     print(_cmpd, _synDrugs, _synPots)
                        #     if _cmpd in _synDrugs:
                        #         _drug, _pot = (_cmpd, None)
                        #     elif _cmpd in _synPots:
                        #         _drug, _pot = (None, _cmpd)
                        #     else:
                        #         _drug, _pot = (None, None)
                        # print(CmpBatchSet, _drug, _pot)
                    #print(SynData)
                #
                # 
                #
                
                # grpData = djTP.wells_df.groupby('cmpbatch_sets') 
                # logger.debug(f" [{djTP.plate_id}] --2-- {djTP.result_type} {djTP.assay_id} {djTP.n_samples} Samples, {djTP.n_inhibitions} Inhibitions")
                # for CmpBatchSet,DRData in grpData:
                #     if CmpBatchSet:
                #         logger.debug(f" [{djTP.plate_id}] --2-- {CmpBatchSet} {len(DRData)}")
                #         # Use CmpBatch w/o SetID
                #         CmpBatch = DRData['cmpbatches'].unique()[0]
                #         logger.debug(f" [{djTP.plate_id}] --3-- {CmpBatch}")
                #         djDR = Synergy().init_data(CmpBatch,DRData[['well_id','conc_lst','inhibition','mscore','conc_unit_lst']],djTP)
                        
                        # _dr.init_data(grpid,grpdf[['well_id','conc_lst','inhibition','conc_unit_lst']],djTP)
                        # if djDR:
                        #     djDR.calc_doseresponse()
                        #     djTP.n_doseresponses += 1
                        #     #print(f" [{djDR.testwell_id}] {str(djDR)}")
                        #     djDR.doseresponse_to_assaydata()
                        #     #print(f" {djDR.assaydata.data_quality} ({djDR.assaydata.data_comment}) ")
                        #     if upload:
                        #         djDR.save_assaydata(overwrite=overwrite)
                        #     #print(repr(djDR))
                        #     DR_Lst.append(djDR)

                if upload:
                    djTP.save()
                
                if verbose>0:
                    logger.info(f" [{djTP.plate_id}] {djTP.result_type} -> {djTP.n_synergy}")
            else:
                logger.info(f" [{djTP.plate_id}] {djTP.result_type} Either no Samples ({djTP.n_samples}) or no Inhibitions ({djTP.n_inhibitions})")

    return(DR_Lst)