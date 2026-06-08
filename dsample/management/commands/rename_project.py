import os
import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dsample.models import Project
from applib.django.foreignkey import rename_ForeignKey

class Command(BaseCommand):
    help = "Rename Projects"

    def add_arguments(self, parser):
        parser.add_argument("--oldid")
        parser.add_argument("--newid")
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--delete",action="store_true",default=False)
        
    def handle(self, *args, **options):
        
        valLog = Validation_Log("Rename_Project")
        
        print(f" [Rename Project] [{options['oldid']}] -> [{options['newid']}] - Upload : {options['upload']} Delete Old: {options['delete']}")
        rename_ForeignKey(Project,options['oldid'],options['newid'],upload=options['upload'],delete=options['delete'])
        
        
# OldID NewID
# PMMV_02 MMV_02
# PMMV_03 MMV_03
# PMMV_04 MMV_04
# PGARDP_01 GARDP_01
# PGARDP_02 GARDP_02
# PHIPS_E01 HIPS_E01
# PHIPS_E02 HIPS_E02
# PHIPS_C01 HIPS_C01
# PHIPS_C02 HIPS_C02
# PENAM_01 ENAMINE_01
# PC_DMSO -> SOLVENT
# PC_001 CTRL_ABX
# PC_002 CTRL_AFX
# PC_003 CTRL_GEN
# PMC_99976 G99976
# PMC_99975 G99975
# PMC_99974 G99974
# PMC_99969 G99969
# PMC_99967 G99967
# PMC_99964 G99964
# PMC_99962 G99962
# PMC_99949 G99949
# PMC_0052 G0052_Protac
# PMC_0051 G0051_FtsZ
# PMC_0050 G0050_GNegHits
# PMC_0048 G0048_Cardlipin
# PMC_0047 G0047_TyrInhib
# PMC_0046 G0046
# PMC_0045 G0045_StrainProfile
# PMC_0043 G0043_Potentiator
# PMC_0037 G0037_Arenicin
# PMC_0036 G0036_Genomics
# PMC_0033 G0033_ClickAB
# PMC_0029 G0029_Anaer
# PMC_0027 G0027_FluoroProbes
# PMC_0002 G0002_Col

    # ABASE_STUDYID = {
    # '010_GPCR':'G0010_GPCR',
    # 'C001_NR4A':'', 'G00_General':'', 'G03_ChemLib':'', 'C008_PSAA':'', 'C002_SOX':'',
    # 'G01_Antibact':'', '026_Antibiotic':'', '011_DSB':'', '003_TB' 'G04_FragLib'
    # '032_NLRP3' '001_Van' 'C014_IMPDH' '033_ClickAB' '800_CO-ADD' '002'
    # '016_Mirabilin' '013_FQHyb' '021_Carb' '005_Ess' '019_Friulimicin'
    # '004_MembBind' 'C009_GHR' 'C004_GLI' 'C012_FIM' '029' '999_Collab'
    # '014_TransGlycInhibit' '027_FluoroProbes' '034_hBD2' '052_TargetDegard'
    # '048' '053_Selenium'        
    # }
