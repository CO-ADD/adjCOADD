import sys, os
import datetime
import numpy as np
import pandas as pd
import scipy.optimize as opt
from decimal import Decimal

from apputil.models import Dictionary
from applib.bio.bio_data import DR_Range, DR2Value
from applib.bio.doseresponse import DoseResponse
from dsample.models import Compound_Batch
from dplate.models import TestPlate, TestWell
from dscreen.models import AssayData_MIC, AssayData_SynergyMIC, AssayData_SynergyIsobol 
from adjcoadd.constants import LIST_SEP, SYN_CLASSES, COMPOUND_SEP, COMPOUNDSET_SEP, AUDITMODEL_FIELDS

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
        
        self.Isobol_lst = []
        self.sumSyn = None
        


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
    @staticmethod
    def _apply_gen_lst(s):
        s['conc_lst_drug'] = [s['conc_drug']]
        s['conc_unit_lst_drug'] = [s['conc_unit_drug']]
        return(s)
    #--------------------------------------------------------------
    
    #--------------------------------------------------------------
    def init_data(self,synDF,DrugID, PotID, synPlate, 
                  RefDRDrug=None, RefDRPot=None, 
                  renameCol= None,
                  CutOff=80, MuSyC=False, Method='min.FICI'):
    #--------------------------------------------------------------
    
        # RENAME_COL = {'conc':'conc_1',     'conc_unit': 'conc_unit_1':, 'COMPOUND_D': 'compound_1': ,
        #               'conc_pot':'conc_2', 'conc_unit_pot': 'conc_unit_2':, 'COMPOUND_P': 'compound_2': ,
        #                'INHIBITION': 'inhibition' }
        
        self.djtp = synPlate
        self.syndf = synDF
        
        if renameCol:
            self.syndf.rename(columns=renameCol, inplace=True )            

        self.drug_id = DrugID
        self.pot_id = PotID
        self.refdr_drug = RefDRDrug
        self.refdr_pot = RefDRPot
        #self.valCol = valCol
        self.cutoff = CutOff
        self.method = Method
        self.musync = MuSyC
        
        self.conc_lst_drug = list(set(self.syndf['conc_drug']))
        self.conc_lst_drug.sort()

        self.conc_lst_pot = list(set(self.syndf['conc_pot']))
        self.conc_lst_pot.sort()
        
        self.syndf = self.syndf.apply(self._apply_gen_lst,axis=1)
                
        #print(f"Checkerboard: {len(self.conc_lst_drug)}x{len(self.conc_lst_pot)}")
        return(self)

    #--------------------------------------------------------------
    def __repr__(self):
        return(f"{self.drug_id} x {self.pot_id}")

    #--------------------------------------------------------------
    def _copy_mic_to_isobol(self,djDR,IsoType=None, SynData=None):
        # Replicate the MIC data to SynergyIsobol for type='Drug' or 'Pot' 
    #--------------------------------------------------------------

        if djDR and IsoType in ['Drug','Pot','SynDrug']:
            # print(f" DEBUG: {[field.name for field in djMIC._meta.get_fields()]}")
            # for f in [field.name for field in djMIC._meta.get_fields()]:
            #     print(f" {f:20s} : {getattr(djMIC,f)}")
            
            djData = djDR.assaydata
            djData.conv_list_to_string() 

            djIsoB = AssayData_SynergyIsobol.get(djData.testplate_id,djData.testwell_id)
            if djIsoB is None:
                djIsoB = AssayData_SynergyIsobol()                        
                djIsoB.testplate_id = djData.testplate_id
                djIsoB.testwell_id = djData.testwell_id
                djIsoB.save_status = 'New'
            else:
                djIsoB.save_status = 'Exists'
                
                
            djIsoB.assay_id = djData.assay_id
            djIsoB.run_id = djData.run_id
                        
            djIsoB.inhibit_max = djData.inhibit_max
            djIsoB.inhibit_min = djData.inhibit_min
            djIsoB.conc_max = djData.conc_max
            djIsoB.conc_min = djData.conc_min
            
            if IsoType == 'Drug':
                djIsoB.cmpbatch_lst_drug = djData.cmpbatch_lst
                djIsoB.cmpbatch_lst_pot = []
                
                djIsoB.synmic_drug = djData.mic
                djIsoB.synmic_unit_drug = djData.mic_unit
                djIsoB.synmic_pot = '-'
                djIsoB.synmic_unit_pot = '-'
                djIsoB.refmic_val_drug = djDR.refmic_val
                djIsoB.refmic_val_pot = -1

                djIsoB.fici_value = 1
                djIsoB.fici_syn = 0
                
            elif IsoType == 'Pot':
                djIsoB.cmpbatch_lst_drug = []
                djIsoB.cmpbatch_lst_pot = djData.cmpbatch_lst
                
                djIsoB.synmic_drug = '-'
                djIsoB.synmic_unit_drug = '-'
                djIsoB.synmic_pot = djData.mic
                djIsoB.synmic_unit_pot = djData.mic_unit
                djIsoB.refmic_val_drug = -1
                djIsoB.refmic_val_pot = djDR.refmic_val

                djIsoB.fici_value = 1
                djIsoB.fici_syn = 0
                
            elif IsoType == 'SynDrug':
                djIsoB.cmpbatch_lst_drug = djData.cmpbatch_lst
                djIsoB.cmpbatch_lst_pot = SynData['pot']
                
                djIsoB.synmic_val_drug = DR2Value(djData.mic)
                djIsoB.synmic_unit_drug = djData.mic_unit
                djIsoB.synmic_val_pot = DR2Value(SynData['conc_pot'])
                djIsoB.synmic_unit_pot = SynData['conc_unit_pot']
                djIsoB.refmic_val_drug = self.refdr_drug.refmic_val
                djIsoB.refmic_val_pot = self.refdr_pot.refmic_val
                
                djIsoB.calc_fici()
            
            djIsoB.synmic = djData.mic
            djIsoB.synmic_unit = djData.mic_unit
            #djIsoB.synmic_prefix = djData.mic_prefix
            
            djIsoB.act_type = djData.act_type
            djIsoB.act_score = djData.act_score
            djIsoB.pscore = djData.pscore

            djIsoB.data_quality = djData.data_quality
            djIsoB.data_comment = djData.data_comment
            djIsoB.valid = djData.valid
            
            djIsoB.analysis = 'pyBioDR'

            #print(repr(djIsoB))
            return(djIsoB)
        
        return(None)
    
    #--------------------------------------------------------------
    def calc_synergy_mic(self):
    #--------------------------------------------------------------
        
        if self.refdr_drug:
            djISO = self._copy_mic_to_isobol(self.refdr_drug,'Drug')
            #print(repr(djISO))
            self.Isobol_lst.append(djISO)
            
        # For each POT conc ---------------------------------------------------------------
        #abCol = {'Conc':valCol['ConcAB'],'Inhib':'INHIBITION','Unit':valCol['UnitAB'],'WellAct':'WELLACT' }

        #print(" [calc_synergy_mic] --------------------------")
        
        grpData = self.syndf.groupby('conc_pot')
        minFICI = None
        minFICI_Value = 99
        minFICI_WellID = 'A00'
        
        for conc_pot,SngData in grpData:
            
            dictSynData = {'conc_pot':conc_pot,
                           'conc_unit_pot':SngData['conc_unit_pot'].unique()[0],
                           'pot':SngData['cmpbatch_pot'].unique().tolist(),
                           'refmic_val_drug':self.refdr_drug.refmic_val,
                           'refmic_val_pot':self.refdr_pot.refmic_val,
                          }
            #print(dictSynData)
            if conc_pot > 0:                                        
                djSynDR = DoseResponse().init_data(self.drug_id,SngData,self.djtp,
                                                  renameCol={'conc_lst_drug':'conc_lst', 'conc_unit_lst_drug':'conc_unit_lst'})
                if djSynDR:
                    djSynDR.calc_doseresponse()
                    djSynDR.doseresponse_to_assaydata()
                    
                    djISO = self._copy_mic_to_isobol(djSynDR,'SynDrug',dictSynData)
                    
                    if djISO.valid > 0:
                        if djISO.synmic[0] == '>':
                            djISO.fici_value = 1
                            djISO.fici_syn = 0
                        
    #                     if '<=' in micAB['MIC']  or '<=' in micPOT['MIC'] :
    #                         IsoB['FICI_SYN'] = -1
    
                        if djISO.fici_value < minFICI_Value:
                            minFICI_Value = djISO.fici_value
                            if minFICI:
                                minFICI_WellID = min(djISO.testwell_id,minFICI.testwell_id)
                            minFICI = djISO
                    else:
                        djISO.fici_value = -9
                        djISO.fici_syn = -9
                    
                    self.Isobol_lst.append(djISO)
                    
        if self.refdr_pot:
            djISO = self._copy_mic_to_isobol(self.refdr_pot,'Pot')
            #print(repr(djISO))
            self.Isobol_lst.append(djISO)

        # minFICI -> sumSyn
        self.sumSyn = AssayData_SynergyMIC.get(minFICI.testplate_id,minFICI.testwell_id)
        if self.sumSyn is None:
            self.sumSyn = AssayData_SynergyMIC()                        
            self.sumSyn.testplate_id = minFICI.testplate_id
            self.sumSyn.testwell_id = minFICI.testwell_id
            self.sumSyn.save_status = 'New'
        else:
            self.sumSyn.save_status = 'Exists'
            
        self.sumSyn.run_id = minFICI.run_id
        self.sumSyn.assay_id = minFICI.assay_id
        
        # Copy values from minFICI        
        for f in [field.name for field in AssayData_SynergyIsobol._meta.get_fields()]:
            if hasattr(self.sumSyn,f) and f not in ['id']+AUDITMODEL_FIELDS:
                setattr(self.sumSyn,f,getattr(minFICI,f))
                
        # Set SynMIC values
        self.sumSyn.checkerboard = f"{len(self.conc_lst_drug)}x{len(self.conc_lst_pot)}"
        self.testwell_id = minFICI_WellID 

    #--------------------------------------------------------------
    def save_isobols(self,overwrite=False):
        # Save self.Isobol_lst
        if len(self.Isobol_lst) > 0:
            for _iso in self.Isobol_lst:
                if _iso.save_status == 'New' or overwrite:
                    _iso.save()

    #--------------------------------------------------------------
    def save_synsum(self,overwrite=False):
        if self.sumSyn:
            if self.sumSyn.save_status == 'New' or overwrite:
                self.sumSyn.save()

    #--------------------------------------------------------------
                        

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
        
        self.n_sumsyn = 0
        self.n_isobol = 0


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
                if verbose > 1:
                    logger.debug(f" [{self.djtp.plate_id}] {self.djtp.result_type} {self.djtp.assay_id} ")
                self.djtp.conv_list_to_string()
                self.djtp.make_wells_df(ListToString=True)

                # Get Synergy Drugs and Potentiators
                # --------------------------------------
                self.syn_drugs = []
                self.syn_pots = []
                self.syn_controls = []
                
                # Get recorded Synergy Drugs and Potentiators from TestPlate
                if self.djtp.synergy_cmpbatches[0]:
                    self.syn_drugs = self.djtp.synergy_cmpbatches[0].split(LIST_SEP)
                if self.djtp.synergy_cmpbatches[1]:
                    self.syn_pots = self.djtp.synergy_cmpbatches[1].split(LIST_SEP)
                if len(self.djtp.synergy_cmpbatches) > 2 and self.djtp.synergy_cmpbatches[2]:
                    self.syn_controls = self.djtp.synergy_cmpbatches[2].split(LIST_SEP)
                else:
                    self.djtp.synergy_cmpbatches.append('')
                    
                #print(f"Synergy Drugs       : {self.syn_drugs}")
                #print(f"Synergy Potentiators: {self.syn_pots}")
                #print(f"Synergy Controls    : {self.syn_controls}")

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
                                              
                #logger.debug(f" [{self.djtp.plate_id}] --1-- Drugs       : {self.syn_drugs} ")
                #logger.debug(f" [{self.djtp.plate_id}] --1-- Potentiators: {self.syn_pots}")
                #logger.debug(f" [{self.djtp.plate_id}] --1-- Controls    : {self.syn_controls}")
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
                dW['well_id'] = r['well_id']
                dW['cmpbatch_drug'] = r['cmpbatch_lst'][drug_id]
                dW['conc_drug'] = r['conc_lst'][drug_id]
                dW['conc_unit_drug'] = r['conc_unit_lst'][drug_id]
                dW['cmpbatch_pot'] = r['cmpbatch_lst'][pot_id]
                dW['conc_pot'] = r['conc_lst'][pot_id]
                dW['conc_unit_pot'] = r['conc_unit_lst'][pot_id]
                #dW['CONC_3'] = conc3
                #dW['CONC_UNIT_3'] = unit3
                #dW['N_COMPOUNDS'] = ncmpd
                dW['inhibition'] = r['inhibition']
                lW.append(dW)
                
                # Check number of Synergy Data
                if dW['conc_drug'] > 0:
                    flSyn[1] += 1
                if dW['conc_pot'] > 0:
                    flSyn[2] += 1
                if dW['conc_drug'] > 0 and dW['conc_pot'] > 0:
                    flSyn[0] += 1

        if 0 not in flSyn:
            mf = pd.DataFrame(lW)
            mf.sort_values(by=['cmpbatch_drug','cmpbatch_pot','conc_drug','conc_pot'], inplace=True)
            return(mf)
        else:
            logger.warning(f"No Synergy Data for [{cDrug} x {cPot}] in Testplate: {self.djtp}")
            return()            

#-------------------------------------------------------------
def get_refmic(djDR):
#--------------------------------------------------------------
    #print("----------------------------")
    ref_val = None
    _multi = 1
    
    # Check if there is a RefMIC from the database
    #-------------------------------------------------
    if djDR.mic_prefix == '>' or djDR.mic_quality in ['Retest','Invalid']:

        _ref_ass = djDR.testplate.assay_id
        _ref_cmp = djDR.testplate.wells[djDR.testwell_id].cmpbatch_lst
        
        _ref_mics = AssayData_MIC.objects.filter(assay_id=_ref_ass,cmpbatch_lst=_ref_cmp,
                                                data_quality='Valid').values_list('mic',flat=True)        
        _ref_range = DR_Range(_ref_mics)
    
        # If there are _ref_mics 
        if _ref_range['nValues'] > 0 :
            _ref_median = _ref_range['Median'].split(COMPOUND_SEP)[0]
            
            # If Median is '<=' try Median with only Values
            if _ref_range['Median'][0] in ['<']:
                _ref_range_valOnly = DR_Range(_ref_mics,ValuesOnly=True)
                if _ref_range_valOnly['nValues'] > 0:
                    _ref_median = _ref_range_valOnly['Median'].split(COMPOUND_SEP)[0]
                    _ref_range = _ref_range_valOnly
                else:
                    _ref_median = _ref_median.replace('<=','')
                
            # If Median is '>' set Muliplier to 2 
            if _ref_median[0] in ['>']:
                _ref_median = _ref_median.replace('>','')
                _multi = 2
            
            # If Median is combination '|' 
            # if COMPOUND_SEP in _ref_median:
            #     _ref_median = _ref_median

            ref_val = _multi*Decimal(_ref_median)

    # Use on Plate MIC
    if ref_val is None:
        _ref_median = djDR.mic.split(COMPOUND_SEP)[0] 
        if _ref_median[0] in ['>']:
            _ref_median = _ref_median.replace('>','')
            _multi = 2
        elif _ref_median[0] in ['<']:
            _ref_median = _ref_median.replace('<=','')

        ref_val = _multi*Decimal(_ref_median)

    return(ref_val)                                 
    

#================================================================================================
def process_testplate_synergy(djTP,upload=False,overwrite=False,verbose=0):
#================================================================================================
    DR_Dict = {}
    ISO_lst = []
    SYN_lst = []
    #djTP = TestPlate.get(PlateID,WellData=True)
    if djTP:
        djSynTP = Synergy_TestPlate().init_data(djTP,verbose=verbose)
        if djSynTP:
            #print('Synergy TestPlate:',djSynTP.djtp.plate_id, djSynTP.djtp.result_type, djSynTP.djtp.assay_id)

            # Single MIC for Drug and Potentiator
            #--------------------------------------
            #print(" ----- Single MIC -------------------------------------")

            grpData = djSynTP.djtp.wells_df.groupby('cmpbatch_sets')            
            for CmpBatchSets,SynData in grpData:
                if CmpBatchSets:
                    CmpBatchs = COMPOUND_SEP.join([TestWell.split_cmpbatch_sets(i)[0] for i in CmpBatchSets.split(COMPOUND_SEP)])
                    if CmpBatchs in djSynTP.syn_drugs+djSynTP.syn_pots+djSynTP.syn_controls:
                        djDR = DoseResponse().init_data(CmpBatchs,SynData[['well_id','conc_lst','inhibition','conc_unit_lst']],djSynTP.djtp)
                        if djDR:
                            djDR.calc_doseresponse()

                            # Get Ref_MIC_Values
                            djDR.refmic_val = get_refmic(djDR)
                                                                                                                    
                            # Prepare AssayData_MIC
                            djDR.doseresponse_to_assaydata() 
                            djSynTP.djtp.n_doseresponses += 1

                            if upload:
                                djDR.save_assaydata(overwrite=overwrite)
                                
                            DR_Dict[CmpBatchs] = djDR
                            

            # Synergy MIC for Drug|Potentiator
            #--------------------------------------
            for _Drug_ID in djSynTP.syn_drugs:
                for _Pot_ID in djSynTP.syn_pots:
                    dfSynDF = djSynTP.get_synergy_data(_Drug_ID,_Pot_ID)
                    
                    if len(dfSynDF) > 0:                    
                        djSyn = Synergy().init_data(dfSynDF,_Drug_ID,_Pot_ID,djSynTP.djtp,
                                                    RefDRDrug=DR_Dict[_Drug_ID],
                                                    RefDRPot=DR_Dict[_Pot_ID], 
                                                    CutOff=80, MuSyC=False, Method='min.FICI')
                        if djSyn:
                            djSyn.calc_synergy_mic() 
                            if djSyn.sumSyn:
                                djSynTP.n_sumsyn += 1
                            if djSyn.Isobol_lst:
                                djSynTP.n_isobol += len(djSyn.Isobol_lst)
                                
                            if upload:
                                djSyn.save_isobols(overwrite=overwrite)
                                djSyn.save_synsum(overwrite=overwrite)
                            
            # Set TestPlate Metadata
            #--------------------------------------
            djTP.synergy_cmpbatches[0] = LIST_SEP.join(djSynTP.syn_drugs)
            djTP.synergy_cmpbatches[1] = LIST_SEP.join(djSynTP.syn_pots)
            djTP.synergy_cmpbatches[2] = LIST_SEP.join(djSynTP.syn_controls)
            djTP.n_synergies = djSynTP.n_sumsyn
            
            if verbose > 1:
                logger.debug(f" [{djSynTP.djtp.plate_id}] {djSynTP.djtp.result_type} {djSynTP.djtp.assay_id} - Syn: {djSynTP.n_sumsyn} - Iso: {djSynTP.n_isobol}")

            if upload:
                djTP.save()

    return(DR_Dict,ISO_lst,SYN_lst)

