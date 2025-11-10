import numpy as np
import pandas as pd
import datetime

#from django_pandas.io import read_frame

from django.db.models import Q

from dsample.models import Project, COADD_Compound
from dplate.models import MasterWell
import logging
logger = logging.getLogger(__name__)


#-----------------------------------------------------------------------------------------
class CmpdPrep_Project():
    """
    Analysis class for Screening data Doseresponse and Single Concentration data
    
    """
    # --------------------------------------------------------------------------------------
    def __init__(self, ProjectID, **kwargs):
    # --------------------------------------------------------------------------------------
        self.n_samples = 0
        self.project_id = ProjectID
        self.project = Project.get(ProjectID)

        self.df_samples = None
        self.dict_samples = {}


        self.DF_COL_CMPD = ['compound_id','compound_code',
                            'reg_mw','reg_amount','reg_amount_unit','reg_solvent',
                            ]

    # --------------------------------------------------------------------------------------
    def get_samples(self):
        self.qryCmpd = COADD_Compound.objects.filter(project_id=self.project
                                                     ).values_list(*self.DF_COL_CMPD)
        self.n_samples = self.qryCmpd.count()
        print(f" [CmpdPrep_Project ] {self.project} {self.n_samples}")
        if self.n_samples > 0:
            self.df_samples = pd.DataFrame(list(self.qryCmpd), columns=self.DF_COL_CMPD)
            self.df_samples = self.df_samples.apply(self.apply_get_barcodes,axis=1)
                
    # --------------------------------------------------------------------------------------
    @staticmethod
    def apply_get_barcodes(s):
        _bc = []
        _mp = []
        _mw = []

        _masterwells = MasterWell.get_barcodes(s['compound_id'])

        if len(_masterwells)>0:
            for key in _masterwells:
                _bc.append(_masterwells[key]['barcode'])
                _mp.append(_masterwells[key]['plate_id'])
                _mw.append(_masterwells[key]['well_id'])
        s['barcode'] = "; ".join(_bc)
        s['masterplate'] = "; ".join(_mp)
        s['masterwell'] = "; ".join(_mw)
        return(s)

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

        if self.n_samples > 0:
            with pd.ExcelWriter(XlFile) as writer:
                if self.n_samples > 0:
                    logger.info(f" [CmpdPrep]     [Samples] {self.df_samples.shape}")
                    self.df_samples.to_excel(writer, sheet_name='Samples')


    # --------------------------------------------------------------------------------------
    def upload_barcodes(self,XlFile=None, **kwargs):
        pass