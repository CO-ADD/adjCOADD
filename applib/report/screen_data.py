import numpy as np
import pandas as pd
import datetime

#from django_pandas.io import read_frame

from django.db.models import Q

from dsummary.models import (Summary_CmpBatch,  Summary_CmpBatch_Doseresp,  Summary_CmpBatch_Inhib,
                             Summary_Structure, Summary_Structure_Doseresp, Summary_Structure_Inhib,)
from dchem.models import Chem_Structure
from dplate.models import TestWell, TestPlate, MasterWell
from dcollab.models import Collab_Group
from dsample.models import Project, COADD_Compound, Library_Compound
from ddrug.models import Drug, VITEK_AST, MIC_COADD
from dscreen.models import AssayData_MIC, AssayData_CC50, AssayData_HC50, Screen_Run, Assay
from applib.bio.bio_data import DR_Range, agg_Inhib, agg_DR, agg_Lst, dr_max_quality, conv_Conc, split_DR, format_DR, DR_GeoMean
from applib.data.dfutils import sort_pivtable_bylevel
from adjcoadd.constants import COMPOUND_SEP

import logging
logger = logging.getLogger(__name__)


#-----------------------------------------------------------------------------------------
class Report_Screening():
    """
    Analysis class for Screening data Doseresponse and Single Concentration data
    
    """
    # --------------------------------------------------------------------------------------
    def __init__(self, **kwargs):
    # --------------------------------------------------------------------------------------

        # - SC Data ------------
        self.COL_TW = [ 'cmpbatch_lst', 'conc_lst','conc_unit_lst','n_cmpbatches',
                        'plate_id__assay_id','inhibition','mscore','act_type','act_score',
                        'plate_id','well_id','plate_id__result_type','plate_id__run_id',
                        ]
        self.DF_COL_SC = [ 'cmpbatch_lst','conc_lst','conc_unit_lst','n_cmpbatches',
                        'assay_id','inhibition','mscore','act_type','act_score',
                        'plate_id','well_id','result_type','run_id',
                        ]

        # - DR Data ------------
        self.COL_MIC  = ['cmpbatch_lst','n_cmpbatches',
                        'testplate_id__assay_id','mic','mic_unit','act_type','act_score','pscore','inhibit_max',                
                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__run_id',
                        'data_quality',
                        ]
        self.COL_CC50 = ['cmpbatch_lst','n_cmpbatches',
                        'testplate_id__assay_id','cc50','cc50_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__run_id',
                        'data_quality',
                        ]
        self.COL_HC50 = ['cmpbatch_lst','n_cmpbatches',
                         'testplate_id__assay_id','hc10','hc50_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','testplate_id__result_type','testplate_id__run_id',
                        'data_quality',
                        ]

        self.DF_COL_DR = ['cmpbatch_lst','n_cmpbatches',
                       'assay_id','dr','dr_unit','act_type','act_score','pscore','inhibit_max',
                        'testplate_id','testwell_id','result_type','run_id',
                        'data_quality',
                        ]

        # - Organisms ---------
        self.ORGANISMS ={ 'COADD' : ['GN_0001','GN_0003','GN_0034','GN_0042','GP_0020','FG_0001','FG_0002'],
                        }
        # - Summary -----------
        self.n_tw = 0
        self.n_mic = 0
        self.n_cc50 = 0
        self.n_hc50 = 0
        self.n_dr = 0
        self.n_sc = 0

        self.n_hcr_sel = 0

        self.n_vitek = 0
        self.n_antibio = 0

        self.n_compounds = 0
        self.n_samples = 0
        self.n_cmpbatch_ids = 0
        self.dict_compounds = {}
        self.dict_samples = {}
        self.list_cmpbatch_ids = []

        self.n_assays = 0
        self.dict_assays = {}

        self.n_testplates = 0
        self.dict_testplates = {}
        self.n_screenruns = 0
        
        self.n_organism_ids = 0
        self.n_cell_ids = 0
        self.list_organism_ids = []
        self.list_cell_ids = []
        
        self.file_name = ''
        
        self.dict_pivtables = {}

        self.verbose = kwargs.get('verbose',0)
        self.valLog = kwargs.get('valLog',None)
    # --------------------------------------------------------------------------------------
    def qry_by_ProjectID(self,ProjectID):
    # --------------------------------------------------------------------------------------
        qryCmpd = COADD_Compound.objects.filter(project_id = ProjectID).values('compound_id','compound_code',)
        self.n_compounds = qryCmpd.count()
        logger.info(f" [Report] ProjectID: {ProjectID} ({self.n_compounds})")
        
        if self.n_compounds > 0:
            self.dict_compounds = {}
            self.list_cmpbatch_ids = []
            for qry in qryCmpd:
                if qry['compound_id'] not in self.dict_compounds:
                    self.dict_compounds[qry['compound_id']] = qry
                    self.dict_compounds[qry['compound_id']]['Source'] = 'COADD'
                    self.list_cmpbatch_ids.append(qry['compound_id'])

            _now = datetime.datetime.now()
            self.file_name = f"Project_{ProjectID}_Summary_{_now:%Y%m%d}"
            
            self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                    cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values_list(*self.COL_MIC)
            self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                    cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values_list(*self.COL_CC50)
            self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                    cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                    testplate_id__plate_quality = 'Valid'                                            
                                    ).values_list(*self.COL_HC50)

            self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                    cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                    plate_id__plate_quality = 'Valid'                                            
                                    ).values_list(*self.COL_TW)
        
        
    # --------------------------------------------------------------------------------------
    def qry_by_RunID(self,RunID_Lst):
    # --------------------------------------------------------------------------------------
        logger.info(f" [Report] RunID: {RunID_Lst} ")
        #print(f" [Report] RunID: {RunID_Lst} ")
        self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                run_id__in = RunID_Lst,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_MIC)
        self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                run_id__in = RunID_Lst,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_CC50)
        self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                run_id__in = RunID_Lst,
                                testplate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_HC50)
        
        self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                plate_id__run_id__in = RunID_Lst,
                                plate_id__plate_quality = 'Valid'                                            
                                ).values_list(*self.COL_TW)
        
        _now = datetime.datetime.now()
        if len(RunID_Lst) == 1:
            self.file_name = self.file_name = f"Run_{RunID_Lst[0]}_Summary_{_now:%Y%m%d}"
        elif len(RunID_Lst) > 1:
            self.file_name = f"Run_{RunID_Lst[0]}_{RunID_Lst[-1]}_Summary_{_now:%Y%m%d}"

    # --------------------------------------------------------------------------------------
    def qry_by_Collaborator(self,CollabGroup,Include_Combination=False):
    # --------------------------------------------------------------------------------------

        if CollabGroup.startswith('CGRP'):
            djGrp = Collab_Group.get(CollabGroup,Code=None, PI_ID=None, Organisation_ID=None)
        else:
            djGrp = Collab_Group.get(None, Code=CollabGroup, PI_ID=None, Organisation_ID=None)

        self.n_compounds = 0
        if djGrp:
            qryCmpd = COADD_Compound.objects.filter(project_id__group_id = djGrp).values('compound_id','compound_code',)
            self.n_compounds = qryCmpd.count()

        logger.info(f" [Report] Collaborator: {CollabGroup} ({self.n_compounds})")


        if self.n_compounds > 0:
            self.dict_compounds = {}
            self.list_cmpbatch_ids = []
            for qry in qryCmpd:
                if qry['compound_id'] not in self.dict_compounds:
                    self.dict_compounds[qry['compound_id']] = qry
                    self.dict_compounds[qry['compound_id']]['Source'] = 'COADD'
                    self.list_cmpbatch_ids.append(qry['compound_id'])
            
            _now = datetime.datetime.now()
            self.file_name = self.file_name = f"Collab_{CollabGroup}_Summary_{_now:%Y%m%d}"

            #print(f" [list_cmpbatch_ids] {len(self.list_cmpbatch_ids)}")
            if Include_Combination:
                # Include any Combinations - SLOW
                self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_MIC)
                self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_CC50)
                self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_HC50)

                self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                        cmpbatch_lst__overlap=self.list_cmpbatch_ids,
                                        plate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_TW)
            else:
                # Only single compounds
                self.qryMIC = AssayData_MIC.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_id__in=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_MIC)
                self.qryCC50 = AssayData_CC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_id__in=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_CC50)
                self.qryHC50 = AssayData_HC50.objects.filter(Q(data_quality = 'Valid') | Q(data_quality = 'Retest'),
                                        cmpbatch_id__in=self.list_cmpbatch_ids,
                                        testplate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_HC50)

                self.qryTW = TestWell.objects.filter(plate_id__result_type='Inhibition', n_cmpbatches__gt = 0,
                                        cmpbatch_id__in=self.list_cmpbatch_ids,
                                        plate_id__plate_quality = 'Valid'                                            
                                        ).values_list(*self.COL_TW)


    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_samples(s):
        s['sample_id'] = COMPOUND_SEP.join(s['cmpbatch_lst'])
        if 'conc_lst' in s:
            s['concs'] = COMPOUND_SEP.join([f"{c}" for c in s['conc_lst']])
        if 'conc_unit_lst' in s:
            s['conc_units'] = COMPOUND_SEP.join(s['conc_unit_lst'])
        return(s)
    
    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_cmpbatch(s):
        s['sample_id'] = COMPOUND_SEP.join(s['cmpbatch_lst'])
        return(s)

    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_dr(s):
        s['dr_max'] = dr_max_quality(s['dr'],s['inhibit_max'],'')
        return(s)
        
    # --------------------------------------------------------------------------------------
    def get_dataframe(self,SC_Only=False, DR_Only=False):
    # --------------------------------------------------------------------------------------
        self.n_samples = 0
        self.n_assays = 0

        self.dict_samples = {}
        self.dict_assays = {}

        _sc = True
        _dr = True
        if DR_Only:
            _sc = False
        elif SC_Only:
            _dr = False

        # - SC Data -------------------------------------------------------
        self.df_sc = None
        if _sc:
            self.n_tw = self.qryTW.count()
            if self.n_tw > 0:
                self.df_sc = pd.DataFrame(list(self.qryTW), columns=self.DF_COL_SC)
                self.df_sc = self.df_sc.apply(self.apply_samples,axis=1)
                logger.info(f" [Report] SC {self.df_sc.shape} [{self.n_tw}] ")
                self.n_sc = self.df_sc.size 

                # - Getting Samples
                for _s in self.df_sc['sample_id'].unique():
                    if _s not in self.dict_samples:
                        self.dict_samples[_s] = {'cmpbatch_lst':_s.split(COMPOUND_SEP)}

                # - Getting Assays
                for _a in self.df_sc['assay_id'].unique():
                    if _a not in self.dict_assays:
                        self.dict_assays[_a] = {'assay_id':_a}

                # - Getting Testplates
                for _p in self.df_sc['plate_id'].unique():
                    if _p not in self.dict_testplates:
                        self.dict_testplates[_p] = {'plate_id':_p}


        # - DR Data -------------------------------------------------------
        dfList = []
        self.df_dr = None
        if _dr:
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
                self.df_dr = self.df_dr.apply(self.apply_samples,axis=1)
                self.df_dr = self.df_dr.apply(self.apply_dr,axis=1)
                self.n_dr = self.df_dr.size
                logger.info(f" [Report] DR: {self.df_dr.shape}  [{self.n_mic} {self.n_cc50} {self.n_hc50}] ")

                # - Getting Samples
                for _s in self.df_dr['sample_id'].unique():
                    if _s not in self.dict_samples:
                        self.dict_samples[_s] = {'cmpbatch_lst':_s.split(COMPOUND_SEP)}

                # - Getting Assays
                for _a in self.df_dr['assay_id'].unique():
                    if _a not in self.dict_assays:
                        self.dict_assays[_a] = {'assay_id':_a}

                # - Getting Testplates
                for _p in self.df_dr['testplate_id'].unique():
                    if _p not in self.dict_testplates:
                        self.dict_testplates[_p] = {'plate_id':_p}

    # --------------------------------------------------------------------------------------
    def get_sample_info(self, Storage_Info=False, Structure_Info=False, Run_Info = False):
    # --------------------------------------------------------------------------------------
        # - Storage Data ------------
        # self.COL_STORAGE    = ['plate_id','well_id', 'barcode','conc_lst','conc_unit_lst']
        # self.DF_COL_STORAGE = ['plate_id','well_id', 'barcode','conc_lst','conc_unit_lst']
        # # - Compounds Data ------------
        # self.COL_CMP = ['cmpbatch_id', 'full_mw','full_mf',
        #                 'batch_source', 'batch_notes']
        # self.DF_COL_CMP = [ 'cmpbatch_id', 'full_mw','full_mf',
        #                 'batch_source', 'batch_notes']

        # self.COL_MCC = ['assay_id','sum_assay_id', 'assay_type',
        #                 'organism_id__organism_name','organism_id__strain_ids','organism_id__strain_code',
        #                 'cell_id__organism_name','cell_id__cell_line',
        #                 ]
        
        _n_samples = [0,0,0]
        _sample_lst = []
        for k in self.dict_samples:
            _dict = {'sample_id':k}

            if k.startswith('MCC_'):
                # MCC sample, check for Drug Info ----------------------------------
                _n_samples[1] += 1
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

            elif k.startswith('C'):
                # CO-ADD  sample ---------------------------------------------------
                _n_samples[0] += 1
                if COADD_Compound.objects.filter(compound_id = k).exists():
                    djCmp = COADD_Compound.get(k)
                    _dict['sample_code'] = djCmp.compound_code
                    _dict['sample_class'] = 'CO-ADD'
                    _dict['project_id'] = str(djCmp.project_id)

                    # Add Storage Info
                    if Storage_Info:
                        qryStorage = MasterWell.objects.filter(cmpbatch_id=k,plate_id__plate_type__in = ['Storage','Master'])
                        n_storage = qryStorage.count()
                        if n_storage > 0:
                            _storage = {'plate_id':[],'well_id':[],'barcode':[],'concs':[],'conc_units':[]}
                            for q in qryStorage:
                                _storage['plate_id'].append(str(q.plate_id))
                                _storage['well_id'].append(q.well_id)

                                if q.conc_lst:
                                    _concs        = COMPOUND_SEP.join([str(x) for x in q.conc_lst if x > 0])
                                    _storage['concs'].append(_concs)
                                else:
                                    _storage['concs'].append('-')

                                if q.conc_unit_lst:
                                    _conc_units   = COMPOUND_SEP.join([str(x) for x in q.conc_unit_lst if x != ""])
                                    _storage['conc_units'].append(_conc_units)
                                else:
                                    _storage['conc_units'].append('-')

                                if q.barcode is None:
                                    _storage['barcode'].append('-')
                                else:
                                    _storage['barcode'].append(q.barcode)
                                
                            _dict['stock_plateid'] = ';'.join(_storage['plate_id'])
                            _dict['stock_wellid'] = ';'.join(_storage['well_id'])
                            _dict['stock_barcode'] = ';'.join(_storage['barcode'])
                            _dict['stock_conc'] = ';'.join(_storage['concs'])
                            _dict['stock_conc_unit'] = ';'.join(_storage['conc_units'])
                        else:
                            _dict['stock_barcode'] = '-'
                            _dict['stock_plateid'] = '-'
                            _dict['stock_wellid'] = '-'
                            _dict['stock_conc'] = 0
                            _dict['stock_conc_unit'] = '-'

                    if Run_Info:
                        #qryTP = TestWell.objects.filter(cmpbatch_id=k)
                        #lstRun_ID = list(TestWell.objects.filter(cmpbatch_id=k).order_by('plate_id__run_id','plate_id__acreated_at').values('plate_id__run_id').distinct())
                        lstRun_ID = list(TestWell.objects.filter(cmpbatch_id=k).order_by('-plate_id__acreated_at').values('plate_id__run_id').distinct())
                        _x = [r for r in lstRun_ID if r['plate_id__run_id'].startswith('H')]
                        if len(_x) >0:
                            _dict['last_runid'] = _x[0]['plate_id__run_id']
                        else:
                            _dict['last_runid'] = '-'

                    # Add Structure Info
                    if Structure_Info:
                        if djCmp.std_smiles != '':
                            _dict['structure_id'] = djCmp.cmpbatch_id.structure_id
                            _dict['smiles'] = djCmp.std_smiles
                        elif djCmp.reg_smiles != '':
                            _dict['structure_id'] = 'NOT REG'
                            _dict['smiles'] = djCmp.reg_smiles
                        else: 
                            _dict['structure_id'] = 'EMPTY'
                            _dict['smiles'] = ''
                                                    
                else:
                    _dict['sample_code'] = "-"
                    _dict['sample_class'] = 'Screen'
                    _dict['project_id'] = '-'
                

            elif k.startswith('LC'):
                # Library sample --------------------------------------------------
                _n_samples[2] += 1

                if Library_Compound.objects.filter(compound_id = k).exists():
                    djCmp = Library_Compound.get(k)
                    _dict['sample_code'] = djCmp.compound_code
                    _dict['sample_class'] = 'Library'
                    _dict['project_id'] = djCmp.library_id
                else:
                    _dict['sample_code'] = "-"
                    _dict['sample_class'] = 'Screen'
                    _dict['project_id'] = '-'


            _sample_lst.append(_dict)
        self.n_samples = len(_sample_lst)
        if  self.n_samples>0:
            self.df_samples = pd.DataFrame(_sample_lst)
            logger.info(f" [Report] Samples: {self.df_samples.shape} {_n_samples} ")
        else:
            logger.warning(f" [Report] No Samples found")
            
    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_assays(s):
        if s['organism_name'] != '-':
            s['assay_org'] = s['organism_name']
        elif s['cell_organism'] != '-':
            s['assay_org'] = s['cell_line']
        else:
            s['assay_org'] = '-'
        return(s)

    # --------------------------------------------------------------------------------------
    def get_assay_info(self):
    # --------------------------------------------------------------------------------------
        # - Assay Data ------------
        self.COL_ASS = ['assay_id','sum_assay_id', 'assay_type','assay_code',
                        'organism_id','organism_id__organism_name','organism_id__strain_ids','organism_id__strain_code',
                        'cell_id','cell_id__organism_name','cell_id__cell_line',
                        ]
        self.DF_COL_ASS = ['assay_id','sum_assay_id', 'assay_type','assay_code',
                        'organism_id','organism_name','strain_ids','strain_code',
                        'cell_id','cell_organism','cell_line',
                        ]
        
        _assay_lst = list(self.dict_assays.keys())
        self.qryAss = Assay.objects.filter(assay_id__in=_assay_lst).values_list(*self.COL_ASS)
        self.n_assays = self.qryAss.count()
        if self.n_assays > 0:
            self.df_assays = pd.DataFrame(list(self.qryAss), columns=self.DF_COL_ASS).fillna('-')
            self.df_assays = self.df_assays.apply(self.apply_assays,axis=1)


            self.list_organism_ids = self.df_assays['organism_id'].unique()
            self.n_organism_ids = len(self.list_organism_ids)
            self.list_cell_ids = self.df_assays['cell_id'].unique()
            self.n_cell_ids = len(self.list_cell_ids)
            logger.info(f" [Report] Assays: {self.n_assays}  [{self.n_organism_ids} {self.n_cell_ids}] ")
        else:
            logger.warning(f" [Report] No Assays found")

    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_testplates(s):
        if 'poscontrol_stats' in s:
            s['posctrl_ave'] = s['poscontrol_stats'][0]
            s['posctrl_std'] = s['poscontrol_stats'][1]

        if 'negcontrol_stats' in s:
            s['negctrl_ave'] = s['negcontrol_stats'][0]
            s['negctrl_std'] = s['negcontrol_stats'][1]

        if 'sample_stats' in s:
            s['sample_ave'] = s['sample_stats'][0]
            s['sample_std'] = s['sample_stats'][1]

        return(s)
            
    # --------------------------------------------------------------------------------------
    def get_testplate_info(self,WithStats=False,WithRunID=True):
    # --------------------------------------------------------------------------------------
        # - Assay Data ------------
        self.COL_TP = ['plate_id','assay_id','run_id','result_type','readout_type',
                       'zfactor','plate_quality','poscontrol_stats','negcontrol_stats','sample_stats',
                       'labware_id__labware_name','labware_id__plate_material','reader',
                        ]
        self.DF_COL_TP = ['plate_id','assay_id','run_id','result_type','readout_type',
                       'zfactor','plate_quality','poscontrol_stats','negcontrol_stats','sample_stats',
                       'labware_name','material','reader',
                        ]

        self.COL_RUN = ['run_id','run_name','run_date','run_conditions','run_type',
                        ]
        self.DF_COL_RUN = ['run_id','run_name','run_date','run_conditions','run_type'
                        ]

        _plate_lst = list(self.dict_testplates.keys())
        self.qryTP = TestPlate.objects.filter(plate_id__in=_plate_lst).values_list(*self.COL_TP)
        self.n_testplates = self.qryTP.count()
        if self.n_testplates > 0:
            self.df_testplates = pd.DataFrame(list(self.qryTP), columns=self.DF_COL_TP).fillna('-')

            if WithRunID:
                _runid_lst = self.df_testplates['run_id'].unique()
                self.qryRun = Screen_Run.objects.filter(run_id__in=_runid_lst).values_list(*self.COL_RUN)
                self.n_screenruns = self.qryRun.count()
                self.df_screenruns = pd.DataFrame(list(self.qryRun), columns=self.DF_COL_RUN).fillna('-')
                logger.info(f" [Report] RunIDs: {self.n_screenruns}  ")

            if WithStats:
                self.df_testplates = self.df_testplates.apply(self.apply_testplates,axis=1)
            logger.info(f" [Report] Testplates: {self.n_testplates}  ")
        else:
            logger.warning(f" [Report] No Testplates found")

    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_vitek(s):
        s['result_type'] = 'Vitek'
        s['run_id'] = s['card_code']
        if '>=' in s['mic']:
            s['mic'] =  s['mic'].replace('>= ','>')
        elif '<=' in s['mic']:
            s['mic'] = s['mic'].replace('<= ','<=')
        s['assay_type'] = s['orgbatch_id'][:7]
        s['assay_org'] =s['organism_name']
        s['dr_max'] = f"{s['mic']} ({s['bp']})"
        return(s)

    # --------------------------------------------------------------------------------------
    def add_vitek_ast(self):
    # --------------------------------------------------------------------------------------
        # - Vitek AST Data ------------
        self.COL_VAST = ['drug_id__drug_name','drug_id__antimicro_class',
                         'card_barcode__orgbatch_id','card_barcode__orgbatch_id__organism_id__organism_name','card_barcode__card_code',
                         'mic','bp_profile',
                        ]
        self.DF_COL_VAST = ['sample_code','sample_class',
                         'orgbatch_id','organism_name','card_code',
                         'mic','bp',
                        ]

        if len(self.list_organism_ids)>0:
            self.qryVAST = (VITEK_AST
                            .objects
                            .filter(card_barcode__orgbatch_id__organism_id__in=self.list_organism_ids)
                            .exclude(mic__exact='')
                            .values_list(*self.COL_VAST)
                            )
            self.n_vitek = self.qryVAST.count()
            if self.n_vitek > 0:
                self.df_vitek = pd.DataFrame(list(self.qryVAST), columns=self.DF_COL_VAST).fillna('-')
                self.df_vitek = self.df_vitek.apply(self.apply_vitek,axis=1)
                logger.info(f" [Report] Vitek AST: {self.df_vitek.shape}  [{self.n_vitek}] ")

    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_antibio(s):
        s['assay_type'] = s['orgbatch_id'][:7]
        s['assay_org'] = s['organism_name']
        return(s)

    @staticmethod
    def apply_antibio_agg(s):
        s['result_type'] = 'BMD'
        s['run_id'] = 'AntiBio'
        # s['assay_type'] = s['orgbatch_id'][:7]
        s['dr_max'] = f"{s['mic']} "
        return(s)

    # --------------------------------------------------------------------------------------
    def add_hcr_selection(self):
    # --------------------------------------------------------------------------------------
        # Generate Selection for HCR only if SC data but no DR data 
        if self.n_sc > 0 and self.n_dr == 0:
            # Filter for 'Active' Samples
            self.n_hcr_sel = 0
            self.df_hcr_sel = None

            _sel_samples = {}
            # Get Sample_ID's for (act_type='A') or (act_type='P' & assay_id like 'GN_')
            for idx,row in self.df_sc.iterrows():
                _sid = row['sample_id']
                if row['act_type'] == 'A':
                    _sel_samples[_sid] = 'A'
                elif row['act_type'] == 'P' and 'GN_' in row['assay_id']:
                    if _sid not in _sel_samples:
                        _sel_samples[_sid] = 'P'
            
            logger.info(f" [SelectHCR] Samples: {len(_sel_samples)}  ")

            # Get pivotSC data for selected samples
            _df_sel_sc = self.df_sc[self.df_sc['sample_id' ].isin([*_sel_samples])]
            _piv_sel_sc = _df_sel_sc.pivot_table(index='sample_id', columns='assay_id', 
                                                        values='act_type',
                                                        aggfunc=lambda x: " ".join(x),)

            # Generate DF for Selected HCR , from df_samples and pivotSC
            _sel_hcr_lst = []
            for _sid in _sel_samples:
                _sample_dict = self.df_samples[self.df_samples['sample_id'] == _sid].to_dict('records')[0]
                _sample_dict['SEL'] = _sel_samples[_sid]
                _sample_piv = _piv_sel_sc.loc[_sid].to_dict()
                _sample_dict.update(_sample_piv)
                _sel_hcr_lst.append(_sample_dict)
                
            self.n_hcr_sel = len(_sel_hcr_lst)
            self.df_hcr_sel = pd.DataFrame(_sel_hcr_lst)

    # --------------------------------------------------------------------------------------
    def add_antibiogram_data(self, RefOrganisms=[]):
    # --------------------------------------------------------------------------------------

        self.COL_ABMIC = ['drug_id__drug_name','drug_id__antimicro_class',
                         'orgbatch_id','orgbatch_id__organism_id__organism_name',
                         'mic','bp_profile',
                        ]
        self.DF_COL_ABMIC = ['sample_code','sample_class',
                         'orgbatch_id','organism_name',
                         'mic','bp',
                        ]
        if len(RefOrganisms)>0:
            _cList = list(self.list_organism_ids) + list(RefOrganisms)
            self.list_organism_ids = list(set(_cList))
            self.n_organism_ids = len(self.list_organism_ids)
        
        if len(self.list_organism_ids)>0:
            self.qryAntiBio = (MIC_COADD
                            .objects
                            .filter(orgbatch_id__organism_id__in=self.list_organism_ids)
                            .exclude(mic__exact='')
                            .values_list(*self.COL_ABMIC)
                            )
            self.n_antibio = self.qryAntiBio.count()
            if self.n_antibio > 0:
                _df_antibio = pd.DataFrame(list(self.qryAntiBio), columns=self.DF_COL_ABMIC).fillna('-')
                _df_antibio = _df_antibio.apply(self.apply_antibio,axis=1)

                showCol = ['sample_code','sample_class','assay_type','assay_org','mic','bp']
                grbyCol = ['sample_code','sample_class','assay_type','assay_org']

                agg_df = (_df_antibio[showCol]
                            .groupby(grbyCol,as_index=False)
                            .agg({'mic':lambda x:DR_Range(x)['Range']})
                            #.aggregate(lambda x: ", ".join(list(np.unique(x)))).sort_values(by=['sample_class'],ascending=True)
                        )

                self.df_antibio = agg_df.apply(self.apply_antibio_agg,axis=1)
                logger.info(f" [Report] Antibiogram : {self.df_antibio.shape}  [{self.n_antibio}] ")


    # --------------------------------------------------------------------------------------
    def gen_pivot_tables(self, PivTables = ['Values','Act'], PivColumns=None, PivRows=None):
    # --------------------------------------------------------------------------------------
        if self.n_sc > 0:
            if not hasattr(self,'df_comb_sc'):
                self.df_comb_sc = self.df_sc.merge(self.df_samples)
            self.df_comb_sc = pd.merge(left=self.df_sc, right=self.df_samples, how= 'left', on='sample_id')
            self.df_comb_sc = pd.merge(left=self.df_comb_sc, right=self.df_assays, how= 'left', on='assay_id')

            self.df_comb_sc = self.df_comb_sc.fillna('-')

            
        if self.n_dr > 0:
            if not hasattr(self,'df_comb_dr'):
                self.df_comb_dr = self.df_dr.merge(self.df_samples)
            self.df_comb_dr = pd.merge(left=self.df_dr, right=self.df_samples, how= 'left', on='sample_id')
            self.df_comb_dr = pd.merge(left=self.df_comb_dr, right=self.df_assays, how= 'left', on='assay_id')

            if self.n_vitek > 0:
                self.df_comb_dr = pd.concat([self.df_comb_dr,self.df_vitek])
            if self.n_antibio > 0:
                self.df_comb_dr = pd.concat([self.df_comb_dr,self.df_antibio])

            self.df_comb_dr = self.df_comb_dr.fillna('-')

        # Setting pivot Rows and Columns
        if PivColumns:
            pivCol = PivColumns
        else:
            pivCol = ['assay_org','assay_type','result_type','run_id']
            if 'AssayID' in PivTables:
                pivCol = ['assay_org','assay_type','assay_id','result_type','run_id']
            pivRow = ['sample_class','sample_code']

        if PivRows:
            pivRow = PivRows
        else:
            pivRow = ['project_id','sample_class','sample_code']


        # -------------------------------------------------------------------------------------------------
        if 'Values' in PivTables:
            if self.n_sc > 0:
                self.piv_sc_ave_inhib = self.df_comb_sc.pivot_table(index=pivRow, 
                                                            columns=pivCol, 
                                                            values='inhibition',
                                                            aggfunc='mean',
                                                            )
                                                            #aggfunc=lambda x: Value_Range(x,floatPrec=1,strSep=";\n")['StrList'])

            if self.n_dr> 0:
                self.piv_dr_drmax = self.df_comb_dr.pivot_table(index=pivRow, 
                                                            columns=pivCol, 
                                                            values='dr_max',
                                                            aggfunc=lambda x: "; ".join(x),
                                                            )
            
            if self.n_dr> 0 and self.n_sc > 0:                                           
                self.piv_values = pd.merge(self.piv_sc_ave_inhib, self.piv_dr_drmax, 'outer', on = pivRow )
                self.dict_pivtables['piv-Values'] = sort_pivtable_bylevel(self.piv_values,0)
            elif self.n_sc > 0 :
                self.dict_pivtables['piv-Values'] = self.piv_sc_ave_inhib
            elif self.n_dr > 0 :
                self.dict_pivtables['piv-Values'] = self.piv_dr_drmax


        # -------------------------------------------------------------------------------------------------
        if 'Act' in PivTables:
            if self.n_sc > 0:
                self.piv_sc_act = self.df_comb_sc.pivot_table(index=pivRow, 
                                                            columns=pivCol, 
                                                            values='act_type',
                                                            aggfunc=lambda x: " ".join(x),
                                                            )
            if self.n_dr > 0:
                self.piv_dr_act = self.df_comb_dr.pivot_table(index=pivRow, 
                                                            columns=pivCol, 
                                                            values='act_type',
                                                            aggfunc=lambda x: " ".join(x),
                                                            )

            if self.n_dr> 0 and self.n_sc > 0:                                           
                self.piv_act = pd.merge(self.piv_sc_act, self.piv_dr_act, 'outer', on = pivRow )
                self.dict_pivtables['piv-Actives'] = sort_pivtable_bylevel(self.piv_act,0) 
            elif self.n_sc > 0 :
                self.dict_pivtables['piv-Actives'] = self.piv_sc_act
            elif self.n_dr > 0 :
                self.dict_pivtables['piv-Actives'] = self.piv_dr_act

 
    # --------------------------------------------------------------------------------------
    def to_datawarrior(self,CsvFile=None, PivColumns=None, PivRows=None):
    # --------------------------------------------------------------------------------------
        # Setting pivot Rows and Columns
        if PivColumns:
            pivCol = PivColumns
        else:
            pivCol = ['assay_code','result_type',]
        pivRow = ['project_id','sample_code','sample_id','smiles']


        if self.n_sc > 0:
            self.datawarrior_sc = self.df_comb_sc.pivot_table(index=pivRow, 
                                                        columns=pivCol, 
                                                        values=['inhibition','mscore'],
                                                        aggfunc={'inhibition':np.mean,
                                                                 'mscore':np.mean}
                                                        )

        if self.n_dr> 0:
            self.datawarrior_dr = self.df_comb_dr.pivot_table(index=pivRow, 
                                                        columns=pivCol, 
                                                        values=['dr_max','pscore'],
                                                        aggfunc={'dr_max':lambda x: " ".join(x), 
                                                                 'pscore':np.mean}
                                                        )
        
        if self.n_dr> 0 and self.n_sc > 0:                                           
            self.datawarrior = pd.merge(self.datawarrior_sc, self.datawarrior_dr, 'outer', on = pivRow )
            self.datawarrior = sort_pivtable_bylevel(self.datawarrior,0)
        elif self.n_sc > 0 :
            self.datawarrior = self.datawarrior_sc
        elif self.n_dr > 0 :
            self.datawarrior = self.datawarrior_dr
            
        self.datawarrior.columns = self.datawarrior.columns.map(' '.join).str.strip()
        self.datawarrior.to_csv(CsvFile)

    # --------------------------------------------------------------------------------------
    def to_excel(self,XlFile=None, Transpose_PivTables=False, verbose=0):
    # --------------------------------------------------------------------------------------
        SHEET_NAME = {
            'piv-Values':'Sum-Values',
            'piv-Actives': 'Sum-ActScore',
            'piv-Data': 'DataWarrior'
        }

        if XlFile is None:
            XlFile = f"{self.file_name}.xlsx" 

        if verbose>0:
            logger.info(f" [Report] Excel --> {XlFile}")

        if self.n_samples > 0:
            with pd.ExcelWriter(XlFile) as writer:
                if self.n_samples > 0:
                    logger.info(f" [Report]     [Samples] {self.df_samples.shape}")
                    self.df_samples.to_excel(writer, sheet_name='Samples')

                if self.n_assays > 0:
                    logger.info(f" [Report]     [Assays] {self.df_assays.shape}")
                    self.df_assays.to_excel(writer, sheet_name='Assays')
                
                if self.n_testplates > 0:
                    COL_EXCLUDE = ['poscontrol_stats','negcontrol_stats','sample_stats']
                    _exp_columns = [c for c in self.df_testplates.columns if c not in COL_EXCLUDE]
                    logger.info(f" [Report]     [TestPlates] {self.df_testplates.shape}")
                    self.df_testplates.to_excel(writer, sheet_name='Testplates',columns=_exp_columns)

                if self.n_screenruns > 0:
                    logger.info(f" [Report]     [Runs] {self.df_screenruns.shape}")
                    self.df_screenruns.to_excel(writer, sheet_name='ScreenRuns')

                if self.n_vitek > 0:
                    logger.info(f" [Report]     [Vitek AST] {self.df_vitek.shape}")
                    self.df_vitek.to_excel(writer, sheet_name='Vitek')

                if self.n_antibio > 0:
                    logger.info(f" [Report]     [AntiBio] {self.df_vitek.shape}")
                    self.df_antibio.to_excel(writer, sheet_name='AntiBio')

                if self.n_hcr_sel > 0:
                    logger.info(f" [Report]     [SelectHCR] {self.df_hcr_sel.shape}")
                    self.df_hcr_sel.to_excel(writer, sheet_name='SelectHCR')

                if self.n_sc > 0:
                    COL_EXCLUDE = ['cmpbatch_lst','conc_lst','conc_unit_lst']
                    _exp_columns = [c for c in self.df_sc.columns if c not in COL_EXCLUDE]
                    _shape = self.df_sc.shape
                    logger.info(f" [Report]     [SC-Data] {_shape}")
                    self.df_sc.to_excel(writer, sheet_name='SC-Data',columns=_exp_columns)

                if self.n_hcr_sel > 0:
                    logger.info(f" [Report]     [HCR Selection] {self.df_hcr_sel.shape}")
                    self.df_hcr_sel.to_excel(writer, sheet_name='HCR Selection')

                if self.n_dr > 0:
                    COL_EXCLUDE = ['cmpbatch_lst']
                    _exp_columns = [c for c in self.df_dr.columns if c not in COL_EXCLUDE]
                    _shape = self.df_dr.shape
                    logger.info(f" [Report]     [DR-Data] {_shape}")
                    self.df_dr.to_excel(writer, sheet_name='DR-Data',columns=_exp_columns)
                
                for k in self.dict_pivtables:
                    logger.info(f" [Report]     [pivTable] {k}")
                    if Transpose_PivTables:
                        self.dict_pivtables[k].T.to_excel(writer, sheet_name=SHEET_NAME[k])
                    else:
                        self.dict_pivtables[k].to_excel(writer, sheet_name=SHEET_NAME[k])
