import sys, os
import datetime
import numpy as np
import pandas as pd
import scipy.optimize as opt

from adjcoadd.constants import COMPOUND_SEP
from applib.bio.bio_data import ActType_DR, pScore, format_DR, dr_max_quality, ActScoreDR_Cutoff
from dsample.models import Compound_Batch
from dscreen.models import AssayData_MIC,AssayData_CC50,AssayData_HC50



import logging
logger = logging.getLogger(__name__)

class DoseResponse():
    
    #--------------------------------------------------------------
    def __init__(self, cutoff=80, 
                 min_dilutions=6,  inhib_limit = 500,
                 inhib_correction = True, inhib_correct_limit = 20,
                 ic50_dmax_cutoff = 40,
                 breakpoint=None, 
                 ):
        self.dr_quality = 'Empty'
        self.dr_type = None
        self.dr_dmax = ''

        self.cutoff = cutoff
        self.min_dilutions = min_dilutions
        self.breakpoint = breakpoint
        self.inhib_limit = inhib_limit
        self.inhib_correct = inhib_correction 
        self.inhib_correct_limit = inhib_correct_limit
        self.ic50_dmax_cutoff = ic50_dmax_cutoff
 
    #--------------------------------------------------------------
    # Initialise DoseRepsonse with
    #   CmpBatches - string of cmpbatches
    #   dfDR - dataframe of dosereponse data
    #   djTestPlate - TestPlate instance
    
    def init_data(self,CmpBatches, dfDR, djTestPlate):

        self.required_columns = ['well_id','conc_lst', 'inhibition','conc_unit_lst']
        # Prepare Data DF
        if set(self.required_columns).issubset(dfDR.columns) :
            
            self.df = dfDR.sort_values(by='conc_lst',ascending=True).reset_index(drop=True)
            self.df = self.df.apply(self.apply_conc,axis=1)
            
            self.dmax = self.df['inhibition'].max(axis=0)
            self.dmin = self.df['inhibition'].min(axis=0)
            self.dave = self.df['inhibition'].mean(axis=0)

            self.cmax = self.df['conc_lst'].max(axis=0)
            self.cmin = self.df['conc_lst'].min(axis=0)
            
            self.n_wells = len(self.df)
            self.n_conc = len(self.df['conc'].unique())
            
            # CmpBatch information - MW of 1st CmpBatch
            self.cmpbatch_lst = CmpBatches.split(COMPOUND_SEP)
            self.n_cmpbatches = len(self.cmpbatch_lst)
            self.cmpbatch_id = Compound_Batch.get(self.cmpbatch_lst[0])
            if self.cmpbatch_id:
                self.full_mw = self.cmpbatch_id.full_mw
            else:
                self.full_mw = 0
            
            
            
            # TestPlate and Well information
            self.testplate = djTestPlate
            self.testwell_id = min(dfDR['well_id'])
            
        else:
            logger.warning(f" [DoseResponse] Missing columns ({dfDR.columns}) [{self.required_columns}]")
            self.df = None
        return(self)

    #--------------------------------------------------------------
    def load_data(self,Plate_ID, Well_ID, Assay_ID):
        pass

    #--------------------------------------------------------------
    def calc_doseresponse(self):
        if 'MIC' in str(self.testplate.result_type): 
            self.MIC()
            self.IC50()
        elif 'CC50' in str(self.testplate.result_type):
            self.IC50() 
        elif 'HC50' in str(self.testplate.result_type):
            self.IC50()
        else:
            logger.warning(f" [DoseResponse] Unknonw ResultType [{self.testplate.result_type}]")
            
    #--------------------------------------------------------------
    def doseresponse_to_assaydata(self,verbose=0):
                
        ASS_FIELDS = { 'MIC': [
                        'cmpbatch_lst','n_cmpbatches','cmpbatch_id',
                        'mic','mic_unit',['mic_skips','skips_active'],
                        ['act_type','mic_act'], ['act_score','mic_act_score'], ['pscore','pmic'],
                        'analysis','n_conc',
                        ['inhibit_max','dmax'], ['inhibit_min','dmin'], ['conc_max','cmax'], ['conc_min','cmin'], 
                        ['data_quality','mic_quality'], ['valid','mic_valid'],
                        'ic50','ic50_unit',['ic50_pscore','pic5'], 
                        'ic50_quality', ['ic50_r2','ic50_fit_r2'], ['ic50_slope','ic50_fit_slope']
                        # ref_mic, ref_mic_chk
                        ],
                    }
        
        ass_key = str(self.testplate.result_type)
        if 'MIC' == ass_key:
            
            self.assaydata_status = 'Exists'
            self.assaydata = AssayData_MIC.get(self.testplate.plate_id,self.testwell_id)
            if self.assaydata is None:
                self.assaydata = AssayData_MIC()
                self.assaydata_status = 'New'
            if verbose > 0:
                logger.info(f" [DoseResponse] {ass_key} ({self.testplate.plate_id}{self.testwell_id}) [{self.assaydata_status}]")
                
        self.assaydata.testplate_id = self.testplate
        self.assaydata.testwell_id = self.testwell_id
        self.assaydata.run_id = self.testplate.run_id
        self.assaydata.assay_id = self.testplate.assay_id
        
        for f in ASS_FIELDS[ass_key]:
            if isinstance(f,list):
                if hasattr(self,f[1]):
                    setattr(self.assaydata,f[0],getattr(self,f[1]))
            else:
                if hasattr(self,f):
                    setattr(self.assaydata,f,getattr(self,f))
                       
        return(self.assaydata)

    #--------------------------------------------------------------
    def save_assaydata(self, overwrite=False):
        if hasattr(self,'assaydata'):
            if self.assaydata:
                self.assaydata.set_defaults_model()
                self.assaydata.validate_fields()
                
                if self.assaydata_status == 'New' or overwrite:    
                    self.assaydata.save()
        
    #--------------------------------------------------------------
    def __str__(self):
        _rstr = []
        if hasattr(self,'mic_dmax'):
            _rstr.append(f"MIC: {self.mic_dmax}")
        if hasattr(self,'ic50_dmax'):
            _rstr.append(f"IC50: {self.ic50_dmax}")
        return("; ".join(_rstr))

    # -- Apply function for Well Activity
    #--------------------------------------------------------------
    @staticmethod
    def apply_conc(s):
        if s['conc_lst']:
            s['conc'] = s['conc_lst'][0]
        else:
            logger.warning(" [DoseResponse] No concentrations ")
        if s['conc_unit_lst']:
            s['conc_unit'] = s['conc_unit_lst'][0]
        else:
            logger.warning(" [DoseResponse] No concentration unit ")
        return(s)
    
    #--------------------------------------------------------------
    @staticmethod
    def apply_active(s, cutoff=80):
        if s['inhibition'] >= cutoff:
            s['active'] = 'A'
        else:
            s['active'] = 'I'
        return(s)

    #-----------------------------------------------------------------------------
    @staticmethod
    def Func_EC50(x, a, b, c, d):
        '''
        Four-parameter log-logistic function
        - a: min response; - b: max response; - c: logEC50; - d: hill slope
        '''
        return(a+(b-a)/(1+np.exp(d*(np.log(x)-np.log(c)))))
    #-----------------------------------------------------------------------------
    @staticmethod
    def Func_IC50(x, c, d):
        '''
        Four-parameter log-logistic function
        - a: 100; - b: 0; - c: logEC50; - d: hill slope
        '''
        a = 100
        b = 0
        return (a+(b-a)/(1+np.exp(d*(np.log(x)-np.log(c)))))
 
    #--------------------------------------------------------------

    #--------------------------------------------------------------

    #===================================================================        
    def MIC(self, cutoff=80):     
    #===================================================================
    
        if self.df is None:
            logger.warning(" [DoseResponse] MIC no data")
            return(None)   
             
        self.df = self.df.apply(self.apply_active,args=(cutoff,),axis=1)

        #----------------------------------------------------------------
        # Evaluate the Well for MIC 
        # - Base on EUCAST - 1 skip -> MIC of lower Conc
        #                  - 2 skips -> reset MIC to higher Conc
        #                  - >2 skips -> Retest

        self.skips_total = 0 # Total #Skips
        self.skips_active = 0 # #Skips within Active Set

        act_Flag = 0
        act_Well= 0
        n_Wells = len(self.df)
        xa1 = 'I'
        xa2 = 'I'
        idx = 1
        for k,v in self.df.iterrows():
            xa0 = v['active']
            if  (xa0 == 'A') and (xa1 == 'I') and (xa2 == 'I') :
                act_Well = idx
                act_Flag = 1
                self.skips_active = 0
            if  (xa0 == 'I') and (act_Flag == 1):
                self.skips_total += 1
                self.skips_active += 1
            xa2 = xa1
            xa1 = xa0
            idx += 1

        # -- Assign correct MIC Value from actWell
        if (act_Well == 0):
            #retMIC = _mic_df.to_dict('records')[nWells-1]
            self.mic_prefix = '>'
            self.mic_well = self.n_wells-1
        elif (act_Well == 1):
            #retMIC = _mic_df.to_dict('records')[act_Well-1]
            self.mic_prefix = '<='
            self.mic_well = act_Well-1
        else:
            #retMIC = _mic_df.to_dict('records')[act_Well-1]
            self.mic_prefix = '='
            self.mic_well = act_Well-1

        self.mic_value = self.df.loc[self.mic_well,'conc_lst']
        self.mic_unit = COMPOUND_SEP.join(self.df.loc[self.mic_well,'conc_unit_lst'])
        self.mic = format_DR(self.mic_prefix,self.mic_value)
        
        self.inhibition_cutoff = cutoff
        #nWells = len(_mic_df)

        self.mic_valid = 1
        self.mic_quality = 'Valid'
        _mic_Comment = []

        # In case n/aSkips  
        if (self.skips_active > 0) and (self.skips_active <= 2):
            self.mic_valid = 1
            self.mic_quality = f'Retest'
            _mic_Comment.append(f"{self.skips_active}/{self.skips_total} Skips")

        elif (self.skips_active > 2):
            self.mic_valid = -1
            self.mic_quality = f"Invalid"
            _mic_Comment.append(f"{self.skips_active} Skips")
        
        # In case CutOff is >80 or <80%  
        if self.inhibition_cutoff < 80:
            self.mic_quality = f'Retest'
            _mic_Comment.append(f"{self.inhibition_cutoff}% Cutoff")
        elif self.inhibition_cutoff > 80:
            _mic_Comment.append(f"{self.inhibition_cutoff}% Cutoff")

        # In case DAve is well outside 0-100 
        if (self.dave > self.inhib_limit) or (self.dave < -self.inhib_limit):
            self.mic_valid = -1
            self.mic_quality = 'Invalid'
            _mic_Comment.append(f"Inhibition")

        # In case not enough dilutions
        if self.n_conc < self.min_dilutions:
            self.mic_valid = 1
            self.mic_quality = 'Retest'
            _mic_Comment.append(f"{self.n_conc} Conc")

        if len(_mic_Comment)>0:
            self.mic_comment = '; '.join(_mic_Comment)
            self.mic_quality += f" ({self.mic_comment})" 
        else:
            self.mic_comment = '-'     
    
        self.analysis = 'pyBioDR'
        
        actID = ActType_DR(self.mic,self.mic_unit,self.dmax)
        self.mic_act_score= ActScoreDR_Cutoff[actID]['Score']
        self.mic_act = ActScoreDR_Cutoff[actID]['Code']

        if self.full_mw > 10:
            self.pmic = pScore(self.mic,self.mic_unit,self.dmax,self.full_mw)
        else:
            self.pmic = -1
        
        self.mic_dmax = dr_max_quality(self.mic,self.dmax,self.mic_quality)
            
        # print(f" [MIC] Inhib: {self.dmin} {self.dmax} Concs: {self.cmin} {self.cmax}")
        # print(f" [MIC] MIC  : {self.mic_dmax} {self.mic_act} {self.mic_act_score} {self.pmic}")  
        # print(f" [MIC]      Skips: {self.skips_active} {self.skips_total}")


    #===================================================================        
    def IC50(self, cutoff=80):     
    #===================================================================

        npConc = np.array(self.df['conc'].to_list(), dtype=float)
        npInhib = np.array(self.df['inhibition'].to_list(), dtype=float)
        
        self.cmax = self.df['conc'].max(axis=0)
        self.cmin = self.df['conc'].min(axis=0)

        # Clip extreme Inhibtion values
        if self.inhib_correct and (self.dmin < -self.inhib_correct-10 or self.dmax > 110+self.inhib_correct):
            npInhib = np.clip(npInhib,-self.inhib_correct,100+self.inhib_correct)
            
        # No activity : DMax < 40%
        if self.dmax < self.ic50_dmax_cutoff:
            self.ic50_prefix = '>'
            self.ic50_value = self.cmax
            self.fit_r2 = 0
            self.ic50_quality = 'Valid'
            self.ic50_valid = 1
        # All Active : DMin > 60%
        elif self.dmin > (100-self.ic50_dmax_cutoff):
            self.ic50_prefix = '<='
            self.ic50_value = self.cmin
            self.fit_r2 = 0
            self.ic50_quality = 'Valid'
            self.ic50_valid = 1
        else:
            try:
                # Func_IC50(x, IC50, slope)
                fit_bounds = ([10**-9, -np.inf],[np.inf, np.inf])
                fit_Coefs, covMat = opt.curve_fit(self.Func_IC50, npConc, npInhib, bounds=fit_bounds)
                res = npInhib - self.Func_IC50(npConc,*fit_Coefs)

                ss_res = np.sum(res**2)
                ss_tot = np.sum((npInhib-np.mean(npInhib))**2)

                self.fit_r2 = round(1 - (ss_res/ss_tot),2)
                self.fit_ic50 = round(fit_Coefs[0],3)
                self.fit_slope = round(fit_Coefs[1],3)
                self.fit_ic10 = round(fit_Coefs[0] * pow(10/90,1/fit_Coefs[1]),3)
                self.fit_ic90 = round(fit_Coefs[0] * pow(90/10,1/fit_Coefs[1]),3)

                if fit_Coefs[0] > self.cmax:
                    self.ic50_prefix = '>'
                    self.ic50_value = self.cmax
                    if self.dmax > self.ic50_dmax_cutoff:
                        self.ic50_prefix = '='
                        self.ic50_value = self.cmax
                        self.fit_r2 = 0
                        
                # Active: fxC50 < Min Tested Concentration
                elif fit_Coefs[0] < self.cmin:
                    self.ic50_prefix = '<='
                    self.ic50_value = self.cmin
                    #retIC50['XC10'] = format_DR('<=',CMin)
                    if self.dmax > self.ic50_dmax_cutoff:
                        self.ic50_prefix = '='
                        self.ic50_value = self.cmin
                        self.fit_r2 = 0
                else:
                    self.ic50_prefix = '='
                    self.ic50_value = self.fit_ic50

                self.ic50_quality = 'Valid'
                self.ic50_valid = 1
                
            # --- No Fit 
            except RuntimeError:
                self.fir_r2 = -1
                if self.dave <= self.ic50_dmax_cutoff:
                    self.ic50_prefix = '>'
                    self.ic50_value = self.cmax
                    self.ic50_quality = 'Retest (NoFit)'
                    self.ic50_valid = 1  
                elif self.dmax > self.ic50_dmax_cutoff:
                    self.ic50_prefix = '='
                    self.ic50_value = self.cmax
                    self.ic50_quality = 'Retest (NoFit)'  
                    self.ic50_valid = 1
                elif self.dave > (100-self.ic50_dmax_cutoff):
                    self.ic50_prefix = '<='
                    self.ic50_value = self.cmin
                    self.ic50_quality = 'Retest (NoFit)'
                    self.ic50_valid = 1
                else:
                    self.ic50_prefix = 'X'
                    self.ic50_value = 0
                    self.ic50_quality = 'Invalid (NoFit)'
                    self.ic50_valid = 0
        
        # If Valid IC50 
        if self.ic50_valid == 1:
            if (self.dave > self.inhib_limit) or (self.dave < -self.inhib_limit):
                self.ic50_quality = 'Invalid (Inhibition)'

            # In case not enough dilutions
            if self.n_conc < self.min_dilutions:
                self.ic50_ = 0
                self.ic50_quality = f"Invalid ({self.n_conc} Conc)"

        # ------------------------------------------------------------------
        # self.ic10 = format_XCFF(retFit['FIT_XC10'],10,self.cmax,self.cmin,self.dmax)
        # self.ic90 = format_XCFF(retFit['FIT_XC90'],90,self.self.cmin,CMin,self.dmax)
        self.ic50_unit = self.df['conc_unit'].unique()[0]

        self.analysis= 'pyBioDR'
        self.ic50 = format_DR(self.ic50_prefix,[self.ic50_value])
        #retIC50['FIT'] = retFit
        
        actID = ActType_DR(self.ic50,self.ic50_unit,self.dmax)
        self.ic50_act_score= ActScoreDR_Cutoff[actID]['Score']
        self.ic50_act = ActScoreDR_Cutoff[actID]['Code']

        if self.full_mw > 10:
            self.pic50 = pScore(self.ic50,self.ic50_unit,self.dmax,self.full_mw)
        else:
            self.pic50 = -1
            
        self.ic50_dmax = dr_max_quality(self.ic50,self.dmax,self.ic50_quality)
        