import os
import pandas as pd
from tqdm import tqdm
from django.core.management.base import BaseCommand, CommandError

from applib.logging.validation_log import Validation_Log
from dsample.models import Project, COADD_Compound
from dchem.models import Chem_Salt

from applib.mol.mol_standardize import SmiStandardizer_DB
from dchem.utils.reg_coadd_compounds import standardize_coadd, checkissues_coadd, regstructure_coadd

#------------------------------------------------
class Command(BaseCommand):
    help = "Register COADD Structures"

    def add_arguments(self, parser):
        parser.add_argument("--01_std",action="store_true",default=False)
        parser.add_argument("--02_chk",action="store_true",default=False)
        parser.add_argument("--03_upd",action="store_true",default=False)
        parser.add_argument("--09_reg",action="store_true",default=False)
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)

    def handle(self, *args, **options):
        
        valLog = Validation_Log("Upload_Project")

        #----------------------------------------
        # reg_smiles -> std_smiles
        #----------------------------------------
        if options["01_std"]:
            qryCmpd = COADD_Compound.objects.all().exclude(std_status='Valid')
            
            nCmpd = qryCmpd.count()
            print(f" [RegCompounds] 01Std [{nCmpd} not Valid] [Upload: {options['upload']} | Overwrite: {options['overwrite']}")
            
            # Declare MolStandardizer with Salt/Ion/Solvent definition from Database
            MolStd = SmiStandardizer_DB(chemdb=Chem_Salt) 
            
            outNumbers = {}            
            for djCmpd in tqdm(qryCmpd.iterator(), total=nCmpd, desc="Processing Compounds"):
                standardize_coadd(djCmpd,MolStd,outNumbers,upload=options['upload'],overwrite=options['overwrite'])
                
            print(f" [RegCompounds] 01Std [{outNumbers}]")

        #----------------------------------------
        # check reg_mw/mf <-> std_mw/mf
        #----------------------------------------
        elif options["02_chk"]:
            qryCmpd = COADD_Compound.objects.all()
            
            nCmpd = qryCmpd.count()
            print(f" [RegCompounds] 02Chk [{nCmpd} All] [Upload: {options['upload']} | Overwrite: {options['overwrite']}")
                        
            outNumbers = {}            
            for djCmpd in tqdm(qryCmpd.iterator(), total=nCmpd, desc="Processing Compounds"):
                checkissues_coadd(djCmpd,outNumbers,upload=options['upload'],overwrite=options['overwrite'])
                
            print(f" [RegCompounds] 02Chk [{outNumbers}]")

        #----------------------------------------
        # update std_mw/mf with changed std_smiles
        #----------------------------------------
        elif options["03_upd"]:
            qryCmpd = COADD_Compound.objects.filter(std_status='Changed') 
            
            nCmpd = qryCmpd.count()
            print(f" [RegCompounds] 03Upd [{nCmpd} Changed] [Upload: {options['upload']} | Overwrite: {options['overwrite']}")
                        
            outNumbers = {}            
            # for djCmpd in tqdm(qryCmpd.iterator(), total=nCmpd, desc="Processing Compounds"):
            #     checkissues_coadd(djCmpd,outNumbers,upload=options['upload'],overwrite=options['overwrite'])
                
            print(f" [RegCompounds] 03Upd [{outNumbers}]")
            
        #----------------------------------------
        # Register Structure from std_smiles
        #----------------------------------------
        elif options["09_reg"]:
            qryCmpd = COADD_Compound.objects.filter(std_status='Valid') 
            
            nCmpd = qryCmpd.count()
            print(f" [RegCompounds] 09Reg [{nCmpd} Valid] [Upload: {options['upload']} | Overwrite: {options['overwrite']}")
                        
            outNumbers = {}            
            for djCmpd in tqdm(qryCmpd.iterator(), total=nCmpd, desc="Processing Compounds"):
                regstructure_coadd(djCmpd,outNumbers,upload=options['upload'],overwrite=options['overwrite'])
                
            print(f" [RegCompounds] 09Reg [{outNumbers}]")
