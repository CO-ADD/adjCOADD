import sys, os
import datetime
import numpy as np
import pandas as pd
from decimal import Decimal, getcontext
from tqdm import tqdm

import logging

from apputil.models import AuditModel
from dchem.models import Chem_Structure
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
#-----------------------------------------------------------------------------

from apputil.models import Dictionary, ApplicationUser


#-----------------------------------------------------------------------------
def chem_structure_export(StructureQry,FileName, FileType='SMI'):
#-----------------------------------------------------------------------------
    
    if FileType=='SMI':
        for q in StructureQry:
            pass
    