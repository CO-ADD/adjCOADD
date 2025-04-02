import numpy as np
import pandas as pd

#from django_pandas.io import read_frame

from django.db.models import Q

from dsummary.models import (Summary_CmpBatch,  Summary_CmpBatch_Doseresp,  Summary_CmpBatch_Inhib,
                             Summary_Structure, Summary_Structure_Doseresp, Summary_Structure_Inhib,)
from dchem.models import Chem_Structure
from dplate.models import TestWell
from dsample.models import Project, COADD_Compound
from ddrug.models import Drug, VITEK_AST
from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run, Assay
from applib.bio.bio_data import DR_Range, conv_Conc, split_DR, format_DR, DR_GeoMean
from adjcoadd.constants import COMPOUND_SEP

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------------------
class Analysis_Screening():
    """
    Analysis class for Screening data Doseresponse and Single Concentration data
    
    """
    # --------------------------------------------------------------------------------------
    def __init__(self):
    # --------------------------------------------------------------------------------------

        # - SC Data ------------
        self.COL_TW = [ 'cmpbatch_lst', 'conc_lst','conc_unit_lst','n_cmpbatches',
                        'plate_id__assay_id','inhibition','mscore','act_type','act_score',
                        'plate_id','well_id','plate_id__result_type',
                        ]
        self.DF_COL_SC = [ 'cmpbatch_lst','conc_lst','conc_unit_lst','n_cmpbatches',
                        'assay_id','inhibition','mscore','act_type','act_score',
                        'plate_id','well_id','result_type',
                        ]

        # - DR Data ------------
        self.COL_MIC  = ['cmpbatch_lst','n_cmpbatches',
                        'testplate_id__assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',                
                        'testplate_id','testwell_id','testplate_id__result_type',
                        ]
        self.COL_CC50 = ['cmpbatch_lst','n_cmpbatches',
                         'testplate_id__assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','testplate_id__result_type',
                        ]
        self.COL_HC50 = ['cmpbatch_lst','n_cmpbatches',
                         'testplate_id__assay_id','hc50','hc50_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','testplate_id__result_type',
                        ]

        self.DF_COL_DR = ['cmpbatch_lst','n_cmpbatches',
                       'assay_id','dr','dr_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','result_type',
                        ]

        
        # - Summary -----------
        self.n_compounds = 0
        self.n_samples = 0
        self.n_assays = 0
        self.n_tw = 0
        self.n_mic = 0
        self.n_cc50 = 0
        self.n_hc50 = 0
        self.n_dr = 0
        self.n_sc = 0

        self.dict_compounds = {}
        self.dict_samples = {}
        self.dict_assays = {}
        self.list_organism_ids = []
        self.list_cmpbatch_ids = []
    # --------------------------------------------------------------------------------------
    def qry_by_ProjectID(self,ProjectID):
    # --------------------------------------------------------------------------------------
        qryCmpd = COADD_Compound.objects.filter(project_id = ProjectID).values('compound_id','compound_code',)
        self.n_compounds = qryCmpd.count()
        logger.info(f" [Analysis] ProjectID: {ProjectID} ({self.n_compounds})")
        
        if self.n_compounds > 0:
            self.dict_compounds = {}
            self.list_cmpbatch_ids = []
            for qry in qryCmpd:
                if qry['compound_id'] not in self.dict_compounds:
                    self.dict_compounds[qry['compound_id']] = qry
                    self.dict_compounds[qry['compound_id']]['Source'] = 'COADD'
                    self.list_cmpbatch_ids.append(qry['compound_id'])

            logger.info(f" [Analysis] ProjectID: {self.n_compounds} ")
            
        self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_MIC)
        self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_CC50)
        self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_HC50)

        self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                plate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_TW)
        
        
    # --------------------------------------------------------------------------------------
    def qry_by_RunID(self,RunID):
    # --------------------------------------------------------------------------------------
        logger.info(f" [Analysis] RunID: {RunID} ")
        self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                run_id = RunID,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_MIC)
        self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                run_id = RunID,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_CC50)
        self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality__contains = 'Retest'),
                                run_id = RunID,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_HC50)
        
        self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                plate_id__run_id = RunID,
                                plate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_TW)

    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_sampleid(s):
        s['sample_id'] = COMPOUND_SEP.join(s['cmpbatch_lst'])
        return(s)
    
    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_vitek(s):
        if '>=' in s['mic']:
            s['mic'] =  s['mic'].replace('>= ','>')
        elif '<=' in s['mic']:
                s['mic'] = s['mic'].replace('<= ','<=')
        return(s)


    # --------------------------------------------------------------------------------------
    def get_dataframe(self):
    # --------------------------------------------------------------------------------------
        self.n_samples = 0
        self.n_assays = 0

        self.dict_samples = {}
        self.dict_assays = {}

        # - SC Data -------------------------------------------------------
        self.df_sc = None
        self.n_tw = self.qryTW.count()
        if self.n_tw > 0:
            self.df_sc = pd.DataFrame(list(self.qryTW), columns=self.DF_COL_SC)
            self.df_sc = self.df_sc.apply(self.apply_sampleid,axis=1)
            logger.info(f" [Analysis] SC {self.df_sc.shape} [{self.n_tw}] ")

            # - Getting Samples
            for _s in self.df_sc['sample_id'].unique():
                if _s not in self.dict_samples:
                    self.dict_samples[_s] = {'cmpbatch_lst':_s.split(COMPOUND_SEP)}

            # - Getting Assays
            for _a in self.df_sc['assay_id'].unique():
                if _a not in self.dict_assays:
                    self.dict_assays[_a] = {'assay_id':_a}


        # - DR Data -------------------------------------------------------
        dfList = []
        self.df_dr = None

        self.n_mic = self.qryMIC.count()
        self.n_cc50 = self.qryCC50.count()
        self.n_hc50 = self.qryHC50.count()
        if self.n_mic > 0:
            dfMIC = pd.DataFrame(list(self.qryMIC), columns=self.DF_COL_DR)
            dfMIC['dr_type'] = 'MIC'
            dfList.append(dfMIC)
        if self.n_cc50 > 0:
            dfCC50 = pd.DataFrame(list(self.qryCC50), columns=self.DF_COL_DR)
            dfCC50['dr_type'] = 'CC50'
            dfList.append(dfCC50)
        if self.n_hc50 > 0:
            dfHC50 = pd.DataFrame(list(self.qryHC50), columns=self.DF_COL_DR)
            dfHC50['dr_type'] = 'HC50'
            dfList.append(dfHC50)
        if dfList:
            self.df_dr = pd.concat(dfList)
            self.df_dr = self.df_dr.apply(self.apply_sampleid,axis=1)
            self.n_dr = self.df_dr.size
            logger.info(f" [Analysis] DR: {self.df_dr.shape}  [{self.n_mic} {self.n_cc50} {self.n_hc50}] ")

            # - Getting Samples
            for _s in self.df_dr['sample_id'].unique():
                if _s not in self.dict_samples:
                    self.dict_samples[_s] = {'cmpbatch_lst':_s.split(COMPOUND_SEP)}

            # - Getting Assays
            for _a in self.df_dr['assay_id'].unique():
                if _a not in self.dict_assays:
                    self.dict_assays[_a] = {'assay_id':_a}

    # --------------------------------------------------------------------------------------
    def get_sample_info(self):
    # --------------------------------------------------------------------------------------
        # - Compounds Data ------------
        # self.COL_CMP = ['assay_id','sum_assay_id', 'assay_type',
        #                 'organism_id__organism_name','organism_id__strain_ids','organism_id__strain_code',
        #                 'cell_id__organism_name','cell_id__cell_line',
        #                 ]
        # self.COL_MCC = ['assay_id','sum_assay_id', 'assay_type',
        #                 'organism_id__organism_name','organism_id__strain_ids','organism_id__strain_code',
        #                 'cell_id__organism_name','cell_id__cell_line',
        #                 ]
        
        # self.DF_COL_CMP = [ 'assay_id','sum_assay_id', 'assay_type',
        #                 'organism_name','strain_ids','strain_code',
        #                 'cell_organism','cell_line',
        #                 ]

        _sample_lst = []
        for k in self.dict_samples:
            _dict = {'sample_id':k}

            # if MCC sample, check for Drug Info
            if 'MCC_' in k:
                _lst = []
                for batch in self.dict_samples[k]['cmpbatch_lst']:
                    _l = batch.split('_')
                    _lst.append("_".join(_l[:2]))
                sample_id = "|".join(_lst)
                self.dict_samples[k]['sample_id'] = sample_id
                if Drug.objects.filter(uq_imb=sample_id).exists():
                    djDrug =   Drug.objects.get(uq_imb=sample_id)               
                    _dict['sample_code'] = djDrug.drug_name
                    _dict['sample_class'] = djDrug.antimicro_class
                    _dict['project_id'] = 'Drugs'
                else:
                    _dict['sample_code'] = "-"
                    _dict['sample_class'] = 'Screen'
                    _dict['project_id'] = '-'

            else:
                pass

            _sample_lst.append(_dict)
        self.n_samples = len(_sample_lst)
        if  self.n_samples>0:
            self.df_samples = pd.DataFrame(_sample_lst)
            

    # --------------------------------------------------------------------------------------
    def get_assay_info(self):
    # --------------------------------------------------------------------------------------
        # - Assay Data ------------
        self.COL_ASS = ['assay_id','sum_assay_id', 'assay_type',
                        'organism_id','organism_id__organism_name','organism_id__strain_ids','organism_id__strain_code',
                        'cell_id__organism_name','cell_id__cell_line',
                        ]
        self.DF_COL_ASS = [ 'assay_id','sum_assay_id', 'assay_type',
                        'organism_id','organism_name','strain_ids','strain_code',
                        'cell_organism','cell_line',
                        ]

        _assay_lst =  self.df_dr['assay_id'].unique()
        self.qrAss = Assay.objects.filter(assay_id__in=_assay_lst).values_list(*self.COL_ASS)
        self.n_assays = self.qrAss.count()
        if self.n_assays > 0:
            self.df_assays = pd.DataFrame(list(self.qrAss), columns=self.DF_COL_ASS).fillna('-')
            self.list_organism_ids = self.df_assays['organism_id'].unique()


    # --------------------------------------------------------------------------------------
    def add_Vitek_AST(self):
    # --------------------------------------------------------------------------------------
        # - Vitek AST Data ------------
        self.COL_VAST = ['drug_id__drug_name','drug_id__drug_codes', 'drug_id__antimicro_class',
                         'card_barcode__orgbatch_id','card_barcode__card_code',
                         'mic','bp_profile',
                        ]
        self.DF_COL_VAST = [ 'drug_name','drug_codes', 'antimicro_class',
                         'orgbatch_id','card_code',
                         'mic','bp_profile',
                        ]

        if len(self.list_organism_ids)>0:
            self.qryVAST = VITEK_AST.objects.filter(card_barcode__orgbatch_id__organism_id__in=self.list_organism_ids).values_list(*self.COL_VAST)
            self.n_vitek = self.qryVAST.count()
            if self.n_vitek > 0:
                self.df_vitek = pd.DataFrame(list(self.qryVAST), columns=self.DF_COL_VAST).fillna('-')
                self.df_vitek = self.df_vitek.apply(self.apply_vitek,axis=1)
                logger.info(f" [Analysis] Vitek AST: {self.df_vitek.shape}  [{self.n_vitek}] ")



    # --------------------------------------------------------------------------------------
    def to_excel(self,XlFile):
    # --------------------------------------------------------------------------------------

        if self.n_dr > 0:
            if not hasattr(self,'df_comb_dr'):
                self.df_comb_dr = self.df_dr.merge(self.df_samples)
            self.df_comb_dr = pd.merge(left=self.df_dr, right=self.df_samples, how= 'left', on='sample_id')
            self.df_comb_dr = pd.merge(left=self.df_comb_dr, right=self.df_assays, how= 'left', on='assay_id')


        drMedian = self.df_comb_dr.pivot_table(index=['sample_class','sample_code'], columns=['organism_name','result_type','assay_id'], values='dr',aggfunc=lambda x: DR_Range(list(x))['Median'])
#        drList = dfMIC_sel.pivot_table(index=['COMPOUND_CODE','COMPOUND'], columns=['ASSAYTYPE_CODE','ASSAY'], values='pDR',aggfunc=lambda x: list(x))
        # if 'L' in Analysis:
        #     print('[sumDR] pivot Data : List DR by [Assay]')
        #     drDList = dfDR_sel.pivot_table(index=['CompoundName','Compound'], columns=['RESULT_TYPE','ASSAYTYPE_CODE','ASSAY'], values='pDR',aggfunc=lambda x: list(x))
        #     drAList = dfDR_sel.pivot_table(index=['CompoundName','Compound'], columns=['RESULT_TYPE','ASSAYTYPE_CODE','ASSAY'], values='ACTIVE',aggfunc=lambda x: list(x))
        # if 'P' in Analysis:
        #     print('[sumDR] pivot Data : DR by [TestPlate]')
        #     drPlateList = dfDR_sel.pivot_table(index=['CompoundName','Compound','TESTWELL_ID'], columns=['RESULT_TYPE','ASSAYTYPE_CODE','ASSAY','TESTPLATE_ID'], values='pDR_LONG',aggfunc=lambda x: x)
        # if 'O' in Analysis:
        #     print('[sumDR] pivot Data : Range DR by [Organism]')
        #     drOrgMedian = dfDR_sel.pivot_table(index=['CompoundCode'], columns=['RESULT_TYPE','ORGANISM'], values='pDR',aggfunc=lambda x: DR_Range(list(x))['Range'])

        # xlFile = os.path.join(OutDir,XlsFileName)
        logger.info(f" [Analysis] Excel : {XlFile}")

        with pd.ExcelWriter(XlFile) as writer:
            # print(f"Datapoints   : {len(dfDR)}") 
            # dfDR.to_excel(writer, sheet_name='DR-Data')

            # print(f"Testplates   : {len(dfTestPlates)}") 
            # dfTestPlates.to_excel(writer, sheet_name='TestPlates')

            # print(f"Compounds    : {len(dfCompounds)}") 
            # dfCompounds.to_excel(writer, sheet_name='Cmpd')

            # print(f"Assays       : {len(dfAssays)}") 
            # dfAssays.to_excel(writer, sheet_name='Assays')
            if self.n_samples > 0:
                logger.info(f" [Analysis] Excel - Samples: {self.df_assays.shape}")
                self.df_samples.to_excel(writer, sheet_name='Samples')

            if self.n_assays > 0:
                logger.info(f" [Analysis] Excel - Assays: {self.df_assays.shape}")
                self.df_assays.to_excel(writer, sheet_name='Assays')
            
            if self.n_vitek > 0:
                logger.info(f" [Analysis] Excel - Vitek AST: {self.df_vitek.shape}")
                self.df_vitek.to_excel(writer, sheet_name='Vitek')


            if self.n_dr > 0:
                _shape = drMedian.shape
                logger.info(f" [Analysis] Excel - DR-Median: {_shape}")
                # if _shape[1]>_shape[0]:
                #     drMedian.T.to_excel(writer, sheet_name='DR-Median')
                # else:
                drMedian.to_excel(writer, sheet_name='DR-Median')

            #micList.to_excel(writer, sheet_name='MIC-List')
            # if 'P' in Analysis:
            #     drPlateList.to_excel(writer, sheet_name='DR-PlateList')
            # if 'O' in Analysis:
            #     drOrgMedian.to_excel(writer, sheet_name='DR-Organism')
            # if 'L' in Analysis:
            #     drDList.to_excel(writer, sheet_name='DR-List')
            #     drAList.to_excel(writer, sheet_name='Active-List')        





