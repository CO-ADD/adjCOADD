from django.db import models
import re

#from django_rdkit import models
from django.db import models
from model_utils import Choices
from sequences import Sequence
from django.core.validators import RegexValidator

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.core.validators import MaxValueValidator, MinValueValidator 
from django.db import transaction, IntegrityError
from django.utils.text import slugify

from apputil.models import AuditModel, Dictionary, ApplicationUser, Document
from applib.data.str_lists import strList_to_List, split_StrList
from applib.bio.bio_data import pScore, ActScore_DR, ActScore_SC
from dsample.models import Compound_Batch

from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------------------------
# Sample/CmpBatch Base/Abstract Models
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
class CmpBatchList_Base(AuditModel):    
#-------------------------------------------------------------------------------------------------

    STRING_FIELDS = ['cmpbatches']
    MAX_CMPBATCHES = 4

    cmpbatches = ""
    cmpbatch_lst = ArrayField(models.CharField(max_length=15, default=""), 
                                 size=MAX_CMPBATCHES, null=True, blank=True, db_index=True, verbose_name = "CmpBatch List")
    n_cmpbatches=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches")

    cmpbatch_id = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id", related_name="%(class)s_cmpbatch_id")

    class Meta:
        abstract = True
        ordering=['cmpbatch_lst']
        # To include in Child Models
        indexes = [
            GinIndex(name="cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="ncmpb_idx", fields=['n_cmpbatches']),
            models.Index(name="cmpbatch_idx", fields=['cmpbatch_lst']),
        ]

    #------------------------------------------------  
    def conv_list_to_string(self):
        if self.cmpbatch_lst:
            _CmpLst = [str(x) for x in self.cmpbatch_lst if x != ""]
            self.cmpbatches   = COMPOUND_SEP.join(_CmpLst)
            self.n_cmpbatches = len(_CmpLst)
        else:
            self.cmpbatches   = ''
            self.n_cmpbatches = 0

    #------------------------------------------------  
    def conv_string_to_list(self):
        self.cmpbatch_lst  = strList_to_List(self.cmpbatches,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def check_cmpbatch_id(self):
        _missing =[]
        for cmpbatch_id in [x for x in self.cmpbatch_lst if x != ""]:
            if not Compound_Batch.exists(cmpbatch_id):
                _missing.append(cmpbatch_id)
        if len(_missing) > 0:
            return({'Error': f"Compound_Batch not found {', '.join(_missing)}"})
        else:
            return(None)

    #------------------------------------------------  
    def set_cmpbatch_id(self,CmpBatchLst=None):
        if CmpBatchLst:
            _CmpLst = [str(x) for x in CmpBatchLst if x != ""]
            setattr(self,'cmpbatch_lst',_CmpLst)
        else:
           _CmpLst = getattr(self,'cmpbatch_lst')
        if _CmpLst: 
            setattr(self,'n_cmpbatches',len(_CmpLst))
            if self.n_cmpbatches == 1:
                self.cmpbatch_id = Compound_Batch.get(_CmpLst[0])
        else:
            setattr(self,'n_cmpbatches',0)
            self.cmpbatch_id = None 

    #------------------------------------------------  
    def __str__(self) -> str:
        return f"{self.cmpbatches}" 

    #------------------------------------------------  
    def __repr__(self) -> str:
        return f"{self.cmpbatches} " 

    #------------------------------------------------
    # Returns an User instance if found by name
    # @classmethod
    # def exists(cls,CmpBatchLst):
    #     if isinstance(CmpBatchLst,str):
    #         CmpBatchLst = [CmpBatchLst]
    #     return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst)).exists()

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def get(cls,CmpBatchLst,verbose=0):
        if isinstance(CmpBatchLst,str):
            CmpBatchLst = [CmpBatchLst]
        _qry = cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst))
        _cnt = _qry.count()
        if _cnt == 1:
            return(_qry[0])
        elif _cnt == 0:
            if verbose:
                logger.warning(f"[CmpBatchList Not Found] {CmpBatchLst} ")
            return(None)
        elif _cnt > 1:
            if verbose:
                logger.warning(f"[CmpBatchList Multiple Exist] {CmpBatchLst} : {_cnt}")
            return(None)
        return(None)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,CmpBatchLst):
        if isinstance(CmpBatchLst,str):
            CmpBatchLst = [CmpBatchLst]
        return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst)).exists()
    
    #------------------------------------------------
    # Clear/Reset cmpbatch data
    def clear_cmpbatch_data(self):
        self.cmpbatches = ""
        self.cmpbatch_lst = []
        self.n_cmpbatches=0
        self.cmpbatch_id = None

#-------------------------------------------------------------------------------------------------
class Sample_Base(CmpBatchList_Base):    
#-------------------------------------------------------------------------------------------------

    STRING_FIELDS = CmpBatchList_Base.STRING_FIELDS + ['concs','conc_units','conc_types'] 

    DICTIONARY_FIELDS = {
        'conc_unit_lst':'Unit_Concentration',
        'conc_type_lst':'Concentration_Type',
    }

    # MAX_CMPBATCHES = 4

    # cmpbatches = ""
    # cmpbatch_lst = ArrayField(models.CharField(max_length=15, default=""), 
    #                              size=MAX_CMPBATCHES, null=True, blank=True, db_index = True, verbose_name = "CmpBatch List")
    # n_cmpbatches=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches")

    concs = ""
    conc_lst = ArrayField(models.DecimalField(max_digits=9, decimal_places=4, default=0), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "Conc List", null=True, blank=True)
    conc_units = ""
    conc_unit_lst = ArrayField(models.CharField(max_length=10, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "ConcUnit List", null=True, blank=True)
    conc_types = ""
    conc_type_lst = ArrayField(models.CharField(max_length=5, default=""), 
                                 size=CmpBatchList_Base.MAX_CMPBATCHES, verbose_name = "ConcType List", null=True, blank=True)

    class Meta:
        abstract = True
        # ordering=['cmpbatch_lst']
        # # To include in Child Models
        # indexes = [
            # models.Index(name="conc_idx", fields=['conc_lst']),
            # models.Index(name="concuni_idx", fields=['conc_unit_lst']),
        # ]

    #------------------------------------------------  
    def conv_list_to_string(self):
        super().conv_list_to_string()
        self.concs        = ''
        self.conc_units   = ''
        self.conc_types   = ''
        
        if self.conc_lst:
            self.concs        = COMPOUND_SEP.join([str(x) for x in self.conc_lst if x > 0])
        if self.conc_unit_lst:
            self.conc_units   = COMPOUND_SEP.join([str(x) for x in self.conc_unit_lst if x != ""])
        if self.conc_type_lst:
            self.conc_types   = COMPOUND_SEP.join([str(x) for x in self.conc_type_lst if x != ""])

    #------------------------------------------------  
    def conv_string_to_list(self):
        super().conv_string_to_list()
        #self.cmpbatch_lst  = strList_to_List(self.cmpbatches,sep=COMPOUND_SEP,size=4,fill="")
        self.conc_lst      = strList_to_List(self.concs,sep=COMPOUND_SEP,size=4,fill=0)
        self.conc_unit_lst = strList_to_List(self.conc_units,sep=COMPOUND_SEP,size=4,fill="")
        self.conc_type_lst = strList_to_List(self.conc_types,sep=COMPOUND_SEP,size=4,fill="")

    #------------------------------------------------  
    def check_conc_unit_dictionary(self):
        _missing = []
        if self.conc_unit_lst:
            for conc_unit in [x for x in self.conc_unit_lst if x != ""]:
                if not Dictionary.exists(self.DICTIONARY_FIELDS['conc_unit_lst'],conc_unit):
                    _missing.append(conc_unit)
        if len(_missing) > 0:
            return({'Error': f"Conc_Unit not found {', '.join(_missing)}"})
        else:
            return(None)

    # #------------------------------------------------  
    # def check_cmpbatch_id(self):
    #     _missing =[]
    #     for cmpbatch_id in [x for x in self.cmpbatch_lst if x != ""]:
    #         if not Compound_Batch.exists(cmpbatch_id):
    #             _missing.append(cmpbatch_id)
    #     if len(_missing) > 0:
    #         return({'Error': f"Compound_Batch not found {', '.join(_missing)}"})
    #     else:
    #         return(None)

    # #------------------------------------------------  
    # def __str__(self) -> str:
    #     return f"{self.cmpbatches}" 

    #------------------------------------------------  
    def __repr__(self) -> str:
        return f"{self.cmpbatches} {self.concs} {self.conc_units}" 

    #------------------------------------------------
    # Clear/Reset cmpbatch data
    def clear_cmpbatch_data(self):
        super().clear_cmpbatch_data()
        self.concs = ""
        self.conc_lst = []
        self.conc_units = ""
        self.conc_unit_lst = []
        self.conc_types = ""
        self.conc_type_lst = []

    # #------------------------------------------------
    # # Returns an User instance if found by name
    # @classmethod
    # def get(cls,CmpBatchLst):
    #     pass

    # #------------------------------------------------
    # # Returns an User instance if found by name
    # @classmethod
    # def exists(cls,CmpBatchLst):
    #     return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst).exists()


#-------------------------------------------------------------------------------------------------
class Comb_CmpBatchList_Base(AuditModel):    
#-------------------------------------------------------------------------------------------------

    COMBINATION_TYPES = ['Drug','Pot']
    STRING_FIELDS = ['cmpbatches_drug','cmpbatches_pot']
    MAX_CMPBATCHES = 2

    # - Drug Drug/Primary 
    cmpbatches_drug = ""
    cmpbatch_lst_drug = ArrayField(models.CharField(max_length=15, default=""), 
                                 size=MAX_CMPBATCHES, null=True, blank=True, db_index=True, verbose_name = "CmpBatch List Drug")
    n_cmpbatches_drug=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches Drug")

    cmpbatch_id_drug = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID Drug", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id_drug", related_name="%(class)s_cmpbatch_id_drug")

    # - Pot Potentiator/Seconadry 
    cmpbatches_pot = ""
    cmpbatch_lst_pot = ArrayField(models.CharField(max_length=15, default=""), 
                                 size=MAX_CMPBATCHES, null=True, blank=True, db_index=True, verbose_name = "CmpBatch List Pot")
    n_cmpbatches_pot=models.SmallIntegerField(default=0, db_index = True,verbose_name = "N CmpBatches Pot")

    cmpbatch_id_pot = models.ForeignKey(Compound_Batch, null=True, blank=True, verbose_name = "CmpBatch ID Pot", on_delete=models.DO_NOTHING,
        db_column="cmpbatch_id_pot", related_name="%(class)s_cmpbatch_id_pot")

    class Meta:
        abstract = True
        ordering=['cmpbatch_lst_drug','cmpbatch_lst_pot']
        # To include in Child Models
        indexes = [
            GinIndex(name="cmp_d_idx",fields=['cmpbatch_lst_drug']),
            models.Index(name="ncmpb_d_idx", fields=['n_cmpbatches_drug']),
            models.Index(name="cmpbatch_d_idx", fields=['cmpbatch_lst_drug']),
            GinIndex(name="cmp_p_idx",fields=['cmpbatch_lst_pot']),
            models.Index(name="ncmpb_p_idx", fields=['n_cmpbatches_pot']),
            models.Index(name="cmpbatch_p_idx", fields=['cmpbatch_lst_pot']),
        ]

    #------------------------------------------------  
    def conv_list_to_string(self):
        if self.cmpbatch_lst_drug:
            _CmpLst = [str(x) for x in self.cmpbatch_lst_drug if x != ""]
            self.cmpbatches_drug   = COMPOUND_SEP.join(_CmpLst)
            self.n_cmpbatches_drug = len(_CmpLst)
        else:
            self.cmpbatches_drug   = ''
            self.n_cmpbatches_drug = 0

        if self.cmpbatch_lst_pot:
            _CmpLst = [str(x) for x in self.cmpbatch_lst_pot if x != ""]
            self.cmpbatches_pot   = COMPOUND_SEP.join(_CmpLst)
            self.n_cmpbatches_pot = len(_CmpLst)
        else:
            self.cmpbatches_pot   = ''
            self.n_cmpbatches_pot = 0

    #------------------------------------------------  
    def conv_string_to_list(self):
        self.cmpbatch_lst_drug  = strList_to_List(self.cmpbatches_drug,sep=COMPOUND_SEP,size=self.MAX_CMPBATCHES,fill="")
        self.cmpbatch_lst_pot   = strList_to_List(self.cmpbatches_pot,sep=COMPOUND_SEP,size=self.MAX_CMPBATCHES,fill="")

    #------------------------------------------------  
    def check_cmpbatch_id(self):
        _missing =[]
        for cmpbatch_id in [x for x in self.cmpbatch_lst_drug if x != ""]:
            if not Compound_Batch.exists(cmpbatch_id):
                _missing.append(cmpbatch_id)
        for cmpbatch_id in [x for x in self.cmpbatch_lst_pot if x != ""]:
            if not Compound_Batch.exists(cmpbatch_id):
                _missing.append(cmpbatch_id)
        if len(_missing) > 0:
            return({'Error': f"Compound_Batch Drug/Pot not found {', '.join(_missing)}"})
        else:
            return(None)

    #------------------------------------------------  
    def set_cmpbatch_id(self,CmpBatchLst=None, type=None):
        
        if type == 'Drug':
            if CmpBatchLst :
                _CmpLst = [str(x) for x in CmpBatchLst if x != ""]
                setattr(self,'cmpbatch_lst_drug',_CmpLst)
            else:
                _CmpLst = getattr(self,'cmpbatch_lst_drug')
                
                if _CmpLst: 
                    setattr(self,'n_cmpbatches_drug',len(_CmpLst))
                    if self.n_cmpbatches == 1:
                        self.cmpbatch_id_drug = Compound_Batch.get(_CmpLst[0])
                else:
                    setattr(self,'n_cmpbatches_drug',0)
                    self.cmpbatch_id_drug = None
        elif type == 'Pot': 
            if CmpBatchLst :
                _CmpLst = [str(x) for x in CmpBatchLst if x != ""]
                setattr(self,'cmpbatch_lst_pot',_CmpLst)
            else:
                _CmpLst = getattr(self,'cmpbatch_lst_pot')
                
                if _CmpLst: 
                    setattr(self,'n_cmpbatches_pot',len(_CmpLst))
                    if self.n_cmpbatches == 1:
                        self.cmpbatch_id_pot = Compound_Batch.get(_CmpLst[0])
                else:
                    setattr(self,'n_cmpbatches_pot',0)
                    self.cmpbatch_id_pot = None 


    #------------------------------------------------  
    def __str__(self) -> str:
        return f"{self.cmpbatches_drug} {COMPOUND_SEP} {self.cmpbatches_pot}" 

    #------------------------------------------------  
    def __repr__(self) -> str:
        return f"{self.cmpbatches} {COMPOUND_SEP} {self.cmpbatches_pot}" 

    #------------------------------------------------
    # Returns an User instance if found by name
    # @classmethod
    # def exists(cls,CmpBatchLst):
    #     if isinstance(CmpBatchLst,str):
    #         CmpBatchLst = [CmpBatchLst]
    #     return cls.objects.filter(cmpbatch_lst__contains=CmpBatchLst, n_cmpbatches = len(CmpBatchLst)).exists()

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def get(cls,CmpBatchLst, type=None, verbose=0):
        if isinstance(CmpBatchLst,str):
            CmpBatchLst = [CmpBatchLst]
        if type == 'Drug':
            _qry = cls.objects.filter(cmpbatch_lst_drug__contains=CmpBatchLst, n_cmpbatches_drug = len(CmpBatchLst))
        elif type == 'Pot':
            _qry = cls.objects.filter(cmpbatch_lst_pot__contains=CmpBatchLst, n_cmpbatches_pot = len(CmpBatchLst))
            
        _cnt = _qry.count()
        if _cnt == 1:
            return(_qry[0])
        elif _cnt == 0:
            if verbose:
                logger.warning(f"[CmpBatchList Not Found] {CmpBatchLst} for {type} ")
            return(None)
        elif _cnt > 1:
            if verbose:
                logger.warning(f"[CmpBatchList Multiple Exist] {CmpBatchLst} for {type} : {_cnt}")
            return(None)
        return(None)

    #------------------------------------------------
    # Returns an User instance if found by name
    @classmethod
    def exists(cls,CmpBatchLst,type=None):
        if type in cls.COMBINATION_TYPES:
            if isinstance(CmpBatchLst,str) :
                CmpBatchLst = [CmpBatchLst]
            if type == 'Drug':
                return cls.objects.filter(cmpbatch_lst_drug__contains=CmpBatchLst, n_cmpbatches_drug = len(CmpBatchLst)).exists()
            elif type == 'Pot':
                return cls.objects.filter(cmpbatch_lst_pot__contains=CmpBatchLst, n_cmpbatches_pot = len(CmpBatchLst)).exists()
        return(None)
    #------------------------------------------------
    # Clear/Reset cmpbatch data
    def clear_cmpbatch_data(self):
        self.cmpbatches_drug = ""
        self.cmpbatch_lst_drug = []
        self.n_cmpbatches_drug=0
        self.cmpbatch_id_drug = None        
        self.cmpbatches_pot = ""
        self.cmpbatch_lst_pot = []
        self.n_cmpbatches_pot=0
        self.cmpbatch_id_pot = None
                