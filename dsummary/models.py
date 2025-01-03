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
from apputil.utils.data import strList_to_List
#from dchem.models import Chem_Structure
from dsample.models import Sample_Base
from adjcoadd.constants import *

import logging
logger = logging.getLogger(__name__)

#-------------------------------------------------------------------------------------------------
# Summary Screening Data Models
#-------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------
class Summary_CmpBatch_Inhib(Sample_Base):
    """
    List of Summary Activity 
    """
#-------------------------------------------------------------------------------------------------
    #from dplate.models import TestPlate

    # Assay Conditions
    assay_id = models.CharField(max_length=25, blank=True, verbose_name = "Assay ID")
    n_assays = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Assay")

    # Activity Summary
    actives = models.CharField(max_length=25, blank=True, verbose_name = "Actives")
    # active_lst
    n_actives = models.SmallIntegerField(default=-1, blank=True, verbose_name = "#Actives")
    act_score_ave = models.DecimalField(default=-1, max_digits=10, decimal_places=2, verbose_name = "Act Score Ave")

    inhibition_ave = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Ave")
    inhibition_std = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Std")
    inhibition_min = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Min")
    inhibition_max = models.DecimalField(default=-1, max_digits=9, decimal_places=3, verbose_name = "Inhibition Max")
    # inhibition_lst
    # inhibition_range = models.CharField(max_length=50, blank=True, verbose_name = "Inhibition Range")
    mscore_ave = models.DecimalField(max_digits=9, decimal_places=3, verbose_name = "MScore Max")

    # conc_lst
    # conc_unit_lst

    # Summary Meta data
    # run_id_lst ArrayField(models.CharField(max_length=15, default="", db_index = True), 
    #                             size=MAX_CMPBATCHES, verbose_name = "CmpBatch List", null=True, blank=True)
    # run_id_date = models.DateField(null=True, blank=True, verbose_name = "Date")

    #------------------------------------------------
    class Meta:
        app_label = 'dsummary'
        db_table = 'sum_cmpbatch_inhib'
        ordering=['assay_id']
        indexes = [
            GinIndex(name="scmpsc_cmp_idx",fields=['cmpbatch_lst']),
            models.Index(name="scmpsc_assid_idx", fields=['assay_id']),
            models.Index(name="scmpsc_nact_idx", fields=['n_actives']),
            models.Index(name="scmpsc_ascr_idx", fields=['act_score_ave']),
            models.Index(name="scmpsc_inhin_idx", fields=['inhibition_ave']),
            models.Index(name="scmpsc_mscr_idx", fields=['mscore_ave']),
        ]
