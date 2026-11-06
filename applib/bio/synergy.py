import sys, os
import datetime
import numpy as np
import pandas as pd
import scipy.optimize as opt

from apputil.models import Dictionary
from applib.bio.bio_data import ActType_DR, pScore, format_DR, dr_max_quality, ActScoreDR_Cutoff
from applib.bio.doseresponse import DoseResponse
from dsample.models import Compound_Batch
from dplate.models import TestPlate, TestWell
from dscreen.models import AssayData_SynergyIsobol, AssayData_SynergySummary
from adjcoadd.constants import SYN_CLASSES, COMPOUND_SEP, COMPOUNDSET_SEP

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

    #--------------------------------------------------------------
    def init_data(self,synDF,cDrug, cPot, micDrug=None, micPot=None,
                valCol={'ConcPOT':'CONC_P',   'UnitPOT':'CONC_UNIT_P',   'CmpdPOT': 'COMPOUND_1',
                        'ConcDRUG' :'CONC_D', 'UnitDRUG' :'CONC_UNIT_D', 'CmpdDRUG' : 'COMPOUND_2',
                        'Inhib':'INHIBITION' },
                CutOff=80, MuSyC=False, Method='min.FICI'):
    #--------------------------------------------------------------
        
        self.syndf = synDF
        self.drug = cDrug
        self.pot = cPot
        self.mic_drug = micDrug
        self.mic_pot = micPot
        self.valCol = valCol
        self.cutoff = CutOff
        self.method = Method
        self.musync = MuSyC

        self.IsoB_Lst = []
        self.sumSyn = {}

    #--------------------------------------------------------------
    def calc_Synergy(self):
    #--------------------------------------------------------------

        IsoB_lst = []

        _concPOT_lst = list(set(synDf[valCol['ConcPOT']]))
        _concPOT_lst.sort()

        minFICI_Value = 999
        minFICI = {}
        # Add AB Data - for Synergy_Plot ---------------------------------------------------------------
        if self.micDrug:
            IsoB = AssayData_SynergyIsobol()
            IsoB['WELL_ID']      = micAB['WELL_ID']

            IsoB['N_CONC_AB']    = micAB['N_CONC']
            IsoB['N_CONC_POT']   = 0
            IsoB['CHECKERBOARD'] = f"{IsoB['N_CONC_AB']}x{IsoB['N_CONC_POT']}"

            IsoB['DMAX']         = micAB['DMAX']
            IsoB['DMIN']         = micAB['DMIN']
            IsoB['CMAX']         = micAB['CMAX']
            IsoB['CMIN']         = micAB['CMIN']
            
            IsoB['COMPOUNDS_AB']    = cAB
            IsoB['N_COMPOUNDS_AB']  = len(IsoB['COMPOUNDS_AB'].split(CmpdSep)) 
            IsoB['COMPOUNDS_POT']   = '-'
            IsoB['N_COMPOUNDS_POT'] = 0
            IsoB['COMPOUNDS']       = IsoB['COMPOUNDS_AB']
            IsoB['N_COMPOUNDS']     = IsoB['N_COMPOUNDS_AB'] + IsoB['N_COMPOUNDS_POT']

            IsoB['SYN_MIC_AB']       = micAB['SYN_MIC']
            IsoB['SYN_UNIT_AB']      = micAB['MIC_UNIT']
            IsoB['SYN_MIC_POT']      = '-'
            IsoB['SYN_UNIT_POT']     = '-'
            
            # IsoB['SYN_MIC']      = zData.merge_lst([IsoB['SYN_AB_MIC'],IsoB['SYN_POT_MIC']],sep=CmpdSep)
            IsoB['SYN_MIC']        = IsoB['SYN_MIC_AB']
            # IsoB['SYN_UNIT']      = zData.merge_lst([IsoB['SYN_UNIT_AB'],IsoB['SYN_UNIT_POT']],sep=CmpdSep)
            IsoB['SYN_UNIT']      = IsoB['SYN_UNIT_AB'] 
            IsoB['SYN_PREFIX'], IsoB['SYN_VALUE'],_ = split_XC50(micAB['MIC'])
            IsoB['SYN_SKIPS']    = micAB['N_SKIPS']
            IsoB['SYN_N_Conc']   = micAB['N_CONC']
            IsoB['DATA_QUALITY']  = 'Valid'
            IsoB['VALID']        = 1
        
            IsoB['PSCORE']       = micAB['PMIC']
            IsoB['ACTSCORE']     = micAB['MIC_ACT_SCORE']
            IsoB['ACTIVE']       = micAB['MIC_ACT']
            IsoB['SYN_IC50']      = micAB['IC50']

            IsoB['FICI_VALUE'] = 1
            IsoB['FICI_SYN'] = 0
            IsoB['ANALYSIS'] = 'pyBioDR'
                    
            self.IsoB_Lst.append(IsoB)
                        
        # For each POT conc ---------------------------------------------------------------
        abCol = {'Conc':valCol['ConcAB'],'Inhib':'INHIBITION','Unit':valCol['UnitAB'],'WellAct':'WELLACT' }

        for concP in _concPOT_lst:
            if concP > 0:
                _synF = synDf[synDf[valCol['ConcPOT']] == concP][synDf[valCol['ConcAB']] != 0]
                #print(_synF)
                
                if (len(_synF)>0):
                    _ab_pot_mic  = calc_MIC_opt(_synF,valCol=abCol)
                    _ab_pot_synmic = get_synMIC(_ab_pot_mic['MIC'],_ab_pot_mic['DMAX'])
                    _ab_pot_ic50 = calc_IC50(_synF,valCol=abCol)
                    _fici, _syn = calc_FICI(micAB['SYN_MIC'],micPOT['SYN_MIC'],_ab_pot_synmic,concP,nCorrect=2)
                    #print(f" {_ab_pot_mic['MIC']} | {concP} - {_ab_pot_mic['DMAX']}")
                    # Collect Isobolagram Data - MIC of AB at different POT conc    
                    IsoB = AssayData_SynergyIsobol()


                    IsoB['WELL_ID']      = _ab_pot_mic['WELL_ID']
                    IsoB['N_CONC_AB']    = len(_synF)
                    IsoB['N_CONC_POT']   = 1
                    IsoB['CHECKERBOARD'] = f"{IsoB['N_CONC_AB']}x{IsoB['N_CONC_POT']}"

                    IsoB['DMAX']         = _ab_pot_mic['DMAX']
                    IsoB['DMIN']         = _ab_pot_mic['DMIN']
                    IsoB['CMAX']         = _ab_pot_mic['CMAX']
                    IsoB['CMIN']         = _ab_pot_mic['CMIN']
                    
                    IsoB['COMPOUNDS_AB']    = cAB
                    IsoB['N_COMPOUNDS_AB']  = len(IsoB['COMPOUNDS_AB'].split(CmpdSep)) 
                    IsoB['COMPOUNDS_POT']   = cPOT
                    IsoB['N_COMPOUNDS_POT'] = len(IsoB['COMPOUNDS_POT'].split(CmpdSep)) 
                    IsoB['COMPOUNDS']       = zData.merge_lst([IsoB['COMPOUNDS_AB'],IsoB['COMPOUNDS_POT']],sep=CmpdSep)
                    IsoB['N_COMPOUNDS']     = IsoB['N_COMPOUNDS_AB'] + IsoB['N_COMPOUNDS_POT']
                    IsoB['SYN_MIC_AB']       = _ab_pot_mic['MIC']
                    IsoB['SYN_UNIT_AB']      = _ab_pot_mic['MIC_UNIT']
                    IsoB['SYN_MIC_POT']      = concP
                    IsoB['SYN_UNIT_POT']     = _ab_pot_mic[valCol['UnitPOT']]
                    
                    IsoB['SYN_MIC']      = zData.merge_lst([IsoB['SYN_MIC_AB'],IsoB['SYN_MIC_POT']],sep=CmpdSep)
                    IsoB['SYN_UNIT']     = zData.merge_lst([IsoB['SYN_UNIT_AB'],IsoB['SYN_UNIT_POT']],sep=CmpdSep)
                    IsoB['SYN_PREFIX'], IsoB['SYN_VALUE'],_ = split_XC50(_ab_pot_mic['MIC'])
                    IsoB['SYN_SKIPS']    = _ab_pot_mic['N_SKIPS']
                    IsoB['SYN_N_Conc']   = _ab_pot_mic['N_CONC']
                    IsoB['DATA_QUALITY']  = _ab_pot_mic['MIC_QUALITY']
                    IsoB['VALID']    = _ab_pot_mic['MIC_VALID']
                
                    IsoB['PSCORE']       = _ab_pot_mic['PMIC']
                    IsoB['ACTSCORE']     = _ab_pot_mic['MIC_ACT_SCORE']
                    IsoB['ACTIVE']       = _ab_pot_mic['MIC_ACT']
                    IsoB['SYN_IC50']      = _ab_pot_ic50['IC50']

                    IsoB['FICI_VALUE'] = 1
                    IsoB['FICI_SYN'] = 0
                    IsoB['ANALYSIS'] = 'pyBioDR'


                    # IsoB['N_AB']         = 
                    
                    # IsoB['SYN_AB']       = _ab_pot_mic[valCol['CmpdAB']]
                    # IsoB['SYN_AB_MIC']   = _ab_pot_mic['MIC']
                    # IsoB['AB_UNIT']      = _ab_pot_mic['MIC_UNIT']
                    # IsoB['SYN_POT']      = _ab_pot_mic[valCol['CmpdPOT']]
                    # IsoB['SYN_POT_MIC']  = cP
                    # IsoB['POT_UNIT']     = _ab_pot_mic[valCol['UnitPOT']]
                    
                    #print(IsoB)
                    #print(IsoB['FICI_VALUE'],minFICI_Value,_ab_pot_mic['MIC_QUALITY'])

                    if _ab_pot_mic['MIC_VALID'] > 0:
                        if IsoB['SYN_PREFIX'] in ['=','<=','<']:
                            IsoB['FICI_SYN'] = _syn
                            IsoB['FICI_VALUE'] = round(_fici,3)
                        else:
                            IsoB['FICI_SYN'] = 0
                            IsoB['FICI_VALUE'] = 1.0
                        
                        if IsoB['FICI_VALUE'] < minFICI_Value:
                            minFICI_Value = IsoB['FICI_VALUE']
                            minFICI = IsoB

                        if '<=' in micAB['MIC']  or '<=' in micPOT['MIC'] :
                            IsoB['FICI_SYN'] = -1
                    else:
                        IsoB['FICI_VALUE'] = -9
                        IsoB['FICI_SYN'] = -9
                        
                    self.IsoB_lst.append(IsoB)

        # Add POT Data - for Synergy_Plot ---------------------------------------------------------------
        if self.micPot:
            IsoB = AssayData_SynergyIsobol()
            IsoB['WELL_ID']      = micPOT['WELL_ID']

            IsoB['N_CONC_AB']    = 0
            IsoB['N_CONC_POT']   = micPOT['N_CONC']
            IsoB['CHECKERBOARD'] = f"{IsoB['N_CONC_AB']}x{IsoB['N_CONC_POT']}"
            IsoB['DMAX']         = micPOT['DMAX']
            IsoB['DMIN']         = micPOT['DMIN']
            IsoB['CMAX']         = micPOT['CMAX']
            IsoB['CMIN']         = micPOT['CMIN']
            
            IsoB['COMPOUNDS_AB']    = '-'
            IsoB['N_COMPOUNDS_AB']  = 0
            IsoB['COMPOUNDS_POT']   = cPOT
            IsoB['N_COMPOUNDS_POT'] = len(IsoB['COMPOUNDS_POT'].split(CmpdSep)) 
            IsoB['COMPOUNDS']       = IsoB['COMPOUNDS_POT']
            IsoB['N_COMPOUNDS']     = IsoB['N_COMPOUNDS_AB'] + IsoB['N_COMPOUNDS_POT']

            IsoB['SYN_MIC_AB']       = '-'
            IsoB['SYN_UNIT_AB']      = '-'
            IsoB['SYN_MIC_POT']      = micPOT['SYN_MIC']
            IsoB['SYN_UNIT_POT']     = micPOT['MIC_UNIT']
            
            # IsoB['SYN_MIC']      = zData.merge_lst([IsoB['SYN_AB_MIC'],IsoB['SYN_POT_MIC']],sep=CmpdSep)
            IsoB['SYN_MIC']        = IsoB['SYN_MIC_POT']
            # IsoB['SYN_UNIT']      = zData.merge_lst([IsoB['SYN_UNIT_AB'],IsoB['SYN_UNIT_POT']],sep=CmpdSep)
            IsoB['SYN_UNIT']      = IsoB['SYN_UNIT_POT'] 
            IsoB['SYN_PREFIX'], IsoB['SYN_VALUE'],_ = split_XC50(micPOT['MIC'])
            IsoB['SYN_SKIPS']    = micPOT['N_SKIPS']
            IsoB['SYN_N_Conc']   = micPOT['N_CONC']
            IsoB['DATA_QUALITY']  = 'Valid'
            IsoB['VALID']    = 1
        
            IsoB['PSCORE']       = micPOT['PMIC']
            IsoB['ACTSCORE']     = micPOT['MIC_ACT_SCORE']
            IsoB['ACTIVE']       = micPOT['MIC_ACT']
            IsoB['SYN_IC50']      = micPOT['IC50']

            IsoB['FICI_VALUE'] = 1
            IsoB['FICI_SYN'] = 0
            IsoB['ANALYSIS'] = 'pyBioDR'

            IsoB_lst.append(IsoB)

                    
        # Summary based on min(FICI) -----------------------------------------------------------
        sumSyn = minFICI.copy()
        sumSyn['ANALYSIS'] = 'pyCastDB'
        sumSyn['WELL_ID'] = min(synDf[synDf[valCol['ConcPOT']] != 0][synDf[valCol['ConcAB']] != 0]['WELL_ID'])
    
        sumSyn['CONC_POT_LST'] =synDf[synDf[valCol['ConcPOT']] != 0][valCol['ConcPOT']].unique()
        sumSyn['CONC_AB_LST']  =synDf[synDf[valCol['ConcAB']] != 0][valCol['ConcAB']].unique()
        sumSyn['N_CONC_POT'] = len(sumSyn['CONC_POT_LST'])
        sumSyn['CHECKERBOARD'] = f"{sumSyn['N_CONC_AB']}x{sumSyn['N_CONC_POT']}"
        
        sumSyn['MIC_AB'] = micAB['MIC']
        sumSyn['MIC_POT'] = micPOT['MIC']

        sumSyn['CONC_AB_MAX'] = synDf[synDf[valCol['ConcAB']]>0][valCol['ConcAB']].max()
        sumSyn['CONC_AB_MIN'] = synDf[synDf[valCol['ConcAB']]>0][valCol['ConcAB']].min()
        sumSyn['CONC_POT_MAX'] = synDf[synDf[valCol['ConcPOT']]>0][valCol['ConcPOT']].max()
        sumSyn['CONC_POT_MIN'] = synDf[synDf[valCol['ConcPOT']]>0][valCol['ConcPOT']].min()

        if sumSyn['N_CONC_AB'] > 6 and sumSyn['N_CONC_POT'] > 6:
            sumSyn['DATA_QUALITY'] = 'Valid'
        elif sumSyn['N_CONC_AB'] > 3 and sumSyn['N_CONC_POT'] > 3:
            sumSyn['DATA_QUALITY'] = 'Retest'
        else:
            sumSyn['DATA_QUALITY'] = f"Invalid ({sumSyn['CHECKERBOARD']})"
            
        if sumSyn['FICI_SYN'] < 0:
            sumSyn['DATA_QUALITY'] = 'Invalid (OutOfRange)'

        # Calculate Syngery based on MuSyc model ------------------------------------------
        if MuSyC:
            if sumSyn['FICI_VALUE'] < 0.6: 
                sMuSyC = calc_Syn_MuSyC(synDf, valCol=valCol, tol=0.5, variant="no_gamma")
                sumSyn['MySyC'] = sMuSyC['MuSyC']
                sumSyn['MySyC_Syn'] = sMuSyC['Synergy']
                sumSyn['MySyC_Ant'] = sMuSyC['Antagonism']
                sumSyn['MySyC_beta'] = sMuSyC['beta']
                sumSyn['MySyC_alpha12'] = sMuSyC['alpha12']
                sumSyn['MySyC_alpha21'] = sMuSyC['alpha21']

        return(sumSyn,IsoB_lst)    



#====================================================================
class Synergy_TestPlate():
#====================================================================
    
    #--------------------------------------------------------------
    def __init__(self, inhibition_cutoff=80, 
                 min_dilutions=6,  inhib_limit = 500, 
                 inhib_correction = True, inhib_correct_limit = 20,
                 ic50_dmax_cutoff = 40,
                 breakpoint=None, 
                 ):

        self.djtp = None
        
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
    def init_data(self,synPlate, verbose=0):
        #
        # Initialise Synergy_TestPlate with TestPlate instance and prepare Synergy Drugs Potentiators Controls
        #
    #--------------------------------------------------------------
        self.djtp = synPlate
        
        if str(self.djtp.result_type) in SYN_CLASSES:
            if self.djtp.n_samples > 0 and self.djtp.n_inhibitions > 0 :
                self.is_syn = True
                self.djtp.n_synergy = 0
                if verbose > 0:
                    logger.info(f" [{self.djtp.plate_id}] {self.djtp.result_type} {self.djtp.assay_id} ")
                self.djtp.conv_list_to_string()
                self.djtp.make_wells_df(ListToString=True)

                # Get Synergy Drugs and Potentiators
                # --------------------------------------
                self.syn_drugs = []
                self.syn_pots = []
                self.syn_controls = []
                
                # Get recorded Synergy Drugs and Potentiators from TestPlate
                if self.djtp.synergy_cmpbatches[0]:
                    self.syn_drugs = self.djtp.synergy_cmpbatches[0].split(';')
                if self.djtp.synergy_cmpbatches[1]:
                    self.syn_pots = self.djtp.synergy_cmpbatches[1].split(';')
                if len(self.djtp.synergy_cmpbatches) > 2 and self.djtp.synergy_cmpbatches[2]:
                    self.syn_controls = self.djtp.synergy_cmpbatches[2].split(';')
                    
                print(f"Synergy Drugs       : {self.syn_drugs}")
                print(f"Synergy Potentiators: {self.syn_pots}")
                print(f"Synergy Controls    : {self.syn_controls}")

                # Complete Synergy Drugs and Potentiators from TestPlate Wells
                #unq_cmpbatches = list(self.djtp.wells_df['cmpbatches'].unique())

                for cmpbatches in self.djtp.wells_df['cmpbatches'].unique():
                    cmpbatch_lst = cmpbatches.split(COMPOUND_SEP)                    
                                        
                    # Get Potentiators based on Drugs
                    for _d in self.syn_drugs:
                        if _d in cmpbatch_lst:
                            cmpbatch_lst.remove(_d)
                            if len(cmpbatch_lst) > 0:
                                _pot = COMPOUND_SEP.join(cmpbatch_lst)
                                if _pot not in self.syn_pots:
                                    self.syn_pots.append(_pot)
                                                                        
                    # Get Drugs based on Potentiators
                    for _p in self.syn_pots:
                        if _p in cmpbatch_lst:
                            cmpbatch_lst.remove(_p)
                            if len(cmpbatch_lst) > 0:
                                _drugs = COMPOUND_SEP.join(cmpbatch_lst)
                                if _drugs not in self.syn_drugs:
                                    self.syn_drugs.append(_drugs)
                                    
                    # Get Controls based on Drugs and Potentiators
                    _druglst = [i for s in [d.split(COMPOUND_SEP) for d in self.syn_drugs] for i in s]
                    _potlst = [i for s in [d.split(COMPOUND_SEP) for d in self.syn_pots] for i in s]
                    _difflst = list(set(cmpbatch_lst) - set(_druglst) - set(_potlst))
                    # print(f"Diff : {_difflst}")
                    for _c in _difflst:
                        if _c and _c not in self.syn_controls:
                            self.syn_controls.append(_c)
                                              
                logger.debug(f" [{self.djtp.plate_id}] --1-- Drugs       : {self.syn_drugs} ")
                logger.debug(f" [{self.djtp.plate_id}] --1-- Potentiators: {self.syn_pots}")
                logger.debug(f" [{self.djtp.plate_id}] --1-- Controls    : {self.syn_controls}")
                return(self)
            
        return(None)

    
    #--------------------------------------------------------------
    def get_synergy_data(self,cDrug,cPot):
        # Extract Synergy Data for Drug|Potentiator from Synergy_TestPlate
    #--------------------------------------------------------------
        lW = []
        flSyn = [0,0,0]
        for i,r in self.djtp.wells_df.iterrows():
            #tmpCmpbatchLst = r['cmpbatch_lst'].copy()
            #tmpCmpbatches = COMPOUND_SEP.join(tmpCmpbatchLst)
            if cDrug in r['cmpbatches'] and cPot in r['cmpbatches']:
                drug_pos = r['cmpbatches'].find(cDrug)
                pot_pos = r['cmpbatches'].find(cPot)
                if drug_pos < pot_pos:
                    drug_id = 0
                    pot_id = 1
                else:
                    drug_id = 1
                    pot_id = 0
                    
                dW = {}
                dW['WELL_ID'] = r['well_id']
                dW['COMPOUND_D'] = r['cmpbatch_lst'][drug_id]
                dW['CONC_D'] = r['conc_lst'][drug_id]
                dW['CONC_UNIT_D'] = r['conc_unit_lst'][drug_id]
                dW['COMPOUND_P'] = r['cmpbatch_lst'][pot_id]
                dW['CONC_P'] = r['conc_lst'][pot_id]
                dW['CONC_UNIT_P'] = r['conc_unit_lst'][pot_id]
                #dW['CONC_3'] = conc3
                #dW['CONC_UNIT_3'] = unit3
                #dW['N_COMPOUNDS'] = ncmpd
                dW['INHIBITION'] = r['inhibition']
                lW.append(dW)
                
                # Check number of Synergy Data
                if dW['CONC_D'] > 0:
                    flSyn[1] += 1
                if dW['CONC_P'] > 0:
                    flSyn[2] += 1
                if dW['CONC_D'] > 0 and dW['CONC_P'] > 0:
                    flSyn[0] += 1

        if 0 not in flSyn:
            mf = pd.DataFrame(lW)
            mf.sort_values(by=['COMPOUND_D','COMPOUND_P','CONC_D','CONC_P'], inplace=True)
            return(mf)
        else:
            return()            

# Process TestPlate for Synergies ------------------------------------------------
# 
def process_testplate_synergy(djTP,upload=False,overwrite=False,verbose=0):
    DR_Lst = []
    #djTP = TestPlate.get(PlateID,WellData=True)
    if djTP:
        djSynTP = Synergy_TestPlate().init_data(djTP,verbose=verbose)
        if djSynTP:
            print('Synergy TestPlate:',djSynTP.djtp.plate_id, djSynTP.djtp.result_type, djSynTP.djtp.assay_id)

            # Single MIC for Drug and Potentiator
            #--------------------------------------
            grpData = djSynTP.djtp.wells_df.groupby('cmpbatch_sets')            
            for CmpBatchSets,SynData in grpData:
                if CmpBatchSets:
                    CmpBatchs = COMPOUND_SEP.join([TestWell.split_cmpbatch_sets(i)[0] for i in CmpBatchSets.split(COMPOUND_SEP)])
                    if CmpBatchs in djSynTP.syn_drugs+djSynTP.syn_pots+djSynTP.syn_controls:
                        djDR = DoseResponse().init_data(CmpBatchs,SynData[['well_id','conc_lst','inhibition','mscore','conc_unit_lst']],djSynTP.djtp)
                        if djDR:
                            djDR.calc_doseresponse()
                            print(repr(djDR))
                            djSynTP.djtp.n_doseresponses += 1
                            
                            djDR.doseresponse_to_assaydata() 
                            #print(f" {djDR.assaydata.data_quality} ({djDR.assaydata.data_comment}) ")
                            if upload:
                                djDR.save_assaydata(overwrite=overwrite)
                            DR_Lst.append(djDR)
                            
            # Synergy MIC for Drug|Potentiator
            #--------------------------------------
            dfSynDF = djSynTP.get_synergy_data(djSynTP.syn_drugs[0],djSynTP.syn_pots[0])
            print(dfSynDF)
            # for pot in djSynTP.syn_pots:
            #     for ab in djSynTP.syn_drugs:
            #         dfSynDF = djSynTP.get_synergy_data(ab,pot)
            #         logger.warning(f" [{djSynTP.djtp.plate_id}] CmpBatchs {CmpBatchs} not in Synergy Drugs, Potentiators or Controls")
        
        # if str(djTP.result_type) in SYN_CLASSES:
        #     if djTP.n_samples > 0 and djTP.n_inhibitions > 0 :
        #         djTP.n_synergy = 0
        #         if verbose > 0:
        #             logger.info(f" [{djTP.plate_id}] {djTP.result_type} {djTP.assay_id} ")
        #         djTP.conv_list_to_string()
        #         djTP.make_wells_df(ListToString=True)

        #         _synDrugs = []
        #         _synPots = []
        #         if djTP.synergy_cmpbatches[0]:
        #             _synDrugs = djTP.synergy_cmpbatches[0].split(';')
        #         if djTP.synergy_cmpbatches[1]:
        #             _synPots = djTP.synergy_cmpbatches[1].split(';')
        #         logger.debug(f" [{djTP.plate_id}] --1-- Drugs: {_synDrugs} - Potentiators: {_synPots}")
                
        #         print(djTP.wells_df)
        #         grpData = djTP.wells_df.groupby('cmpbatch_sets')
                
        #         for CmpBatchSet,SynData in grpData:
        #             if CmpBatchSet:
        #                 # Single MIC for Drug and Potentiator
        #                 #--------------------------------------
        #                 _cmpdsets = CmpBatchSet.split(COMPOUND_SEP)
        #                 if len(_cmpdsets) == 1:
        #                     _cmpbatch = TestWell.split_cmpbatch_sets(_cmpdsets[0])[0]
        #                     djDR = DoseResponse().init_data(_cmpbatch,SynData[['well_id','conc_lst','inhibition','mscore','conc_unit_lst']],djTP)
        #                     if djDR:
        #                         djDR.calc_doseresponse()
        #                         print(repr(djDR))
        #                         djTP.n_doseresponses += 1
                                
        #                         djDR.doseresponse_to_assaydata() 
        #                         #print(f" {djDR.assaydata.data_quality} ({djDR.assaydata.data_comment}) ")
        #                         if upload:
        #                             djDR.save_assaydata(overwrite=overwrite)
        #                         DR_Lst.append(djDR)

        #         for CmpBatchSet,SynData in grpData:
        #             if CmpBatchSet:
        #                 # Synergy MIC for Drug|Potentiator
        #                 #--------------------------------------
        #                 _cmpdsets = CmpBatchSet.split(COMPOUND_SEP)
        #                 if len(_cmpdsets) > 1:
        #                     print(_cmpdsets, _synDrugs, _synPots)
                            
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

            #     if upload:
            #         djTP.save()
                
            #     if verbose>0:
            #         logger.info(f" [{djTP.plate_id}] {djTP.result_type} -> {djTP.n_synergy}")
            # else:
            #     logger.info(f" [{djTP.plate_id}] {djTP.result_type} Either no Samples ({djTP.n_samples}) or no Inhibitions ({djTP.n_inhibitions})")

    return(DR_Lst)