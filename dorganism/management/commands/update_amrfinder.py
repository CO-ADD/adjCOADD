import os, sys
import pandas as pd
import logging
from tqdm import tqdm
from decimal import Decimal

from adjcoadd import settings
from django.core.management.base import BaseCommand, CommandError
from applib.logging.validation_log import Validation_Log

from dorganism.models import Organism
from dgene.models import AMR_Genotype,ID_Sequence
# from dscreen.models import Screen_Run
# from adjcoadd.constants import DR_CLASSES, SYN_CLASSES
# from applib.plate.testplates import add_mother_to_testplate
# from applib.bio.doseresponse import calc_testplate_doseresponse, reset_testplate_assaydata
# from applib.bio.synergy import calc_testplate_synergy


logger = logging.getLogger(__name__)

#================================================================================================
class Command(BaseCommand):
#================================================================================================
    
    command_process = "Update AMR Finder"

    def add_arguments(self, parser):

        parser.add_argument("--orgid",action="store",default=None,help='OrgID')
        parser.add_argument("--all",action="store_true",default=False,help='All')
        
        parser.add_argument("--upload",action="store_true",default=False)
        parser.add_argument("--overwrite",action="store_true",default=False)
        parser.add_argument("--debug",action="store_true",default=False)


    @staticmethod
    def _split_listproperties(ListProp):
        PropDict = {}
        _l = ListProp.split(';')
        for s in _l:
            start = s.find("[")
            end = s.find("]", start + 1)
            if start != -1 and end != -1:
                key = s[start + 1 : end]
                val = s[end+1:]
                PropDict[key] = val.strip()
        return(PropDict)

    @staticmethod
    def _gen_listproperties(PropDict):
        _list =[]
        for key in PropDict:
            _list.append(f"[{key}] {PropDict[key]}")
        return('; '.join(_list))
    
    
    #------------------------------------------------------------------  
    def _process_organism(self,OrgID, *args, **kwargs):
    #------------------------------------------------------------------  
        upload = kwargs.get('upload', False)
        overwrite = kwargs.get('overwrite', False)

        verbose = kwargs.get('verbose', 0)
        valLog = kwargs.get('valLog', None)
        
        
        _Changed = False
        
        djOrg = Organism.get(OrgID)
        
        if djOrg is None:
            logger.error(f" [{self.command_process}] Organism ID {OrgID} not found")
            return()


        # Gen_Property
        #----------------------------------------------------------------------------
        _Valid = True
        if djOrg.gen_property == '':
            propdict = {}
        else:
            propdict = self._split_listproperties(djOrg.gen_property)
            if len(propdict) == 0:
                _Valid = False
                _Changed = False
                logger.error(f" {OrgID} Wrong GenProperty Format [{djOrg.gen_property}] ")
                
        if _Valid:
            qryAMR = AMR_Genotype.objects.filter(seq_id__orgbatch_id__organism_id = djOrg, 
                                                    seq_id__orgbatch_id__qc_status = 'Pass',
                                                    gene_id__gene_type = 'Resistance', 
                                                    ).values_list('amr_method','gene_id__gene_code','gene_id__gene_subtype')
            if qryAMR.count() > 0 :
                amrgene_df = pd.DataFrame(list(qryAMR), columns=['amr_method','gene_code','gene_subtype'])
                if len(amrgene_df)>0:
                    amrfinder_gene_lst = amrgene_df[amrgene_df['amr_method'] == 'AMR Finder']['gene_code'].to_list()                    
                    propdict['AMR Finder'] =  ", ".join(amrfinder_gene_lst)
                    djOrg.gen_property = self._gen_listproperties(propdict)
                    _Changed = True

        # Sero_Clone
        #----------------------------------------------------------------------------
        _Valid = True
        if djOrg.sero_clone == '':
            serodict = {}
        else:
            serodict = self._split_listproperties(djOrg.sero_clone)
            if len(serodict) == 0:
                _Valid = False
                _Changed = False
                logger.error(f" {OrgID} Wrong SeroType Format [{djOrg.sero_clone}] ")
                
        if _Valid:
            qryMLST = ID_Sequence.objects.filter(seq_id__orgbatch_id__organism_id = djOrg, 
                                                    seq_id__orgbatch_id__qc_status = 'Pass',
                                                    seq_filetype = 'FastA', 
                                                    ).values_list('mlst_scheme','mlst_seqtype','mlst_alleles')
            if qryMLST.count() > 0 :
                mlst_df = pd.DataFrame(list(qryMLST), columns=['mlst_scheme','mlst_seqtype','mlst_alleles'])
                mlst_df = mlst_df[mlst_df['mlst_seqtype'] != '-']
                print(mlst_df)
                #amrfinder_gene_lst = amrgene_df[amrgene_df['amr_method'] == 'AMR Finder']['gene_code'].to_list()                    
            
        
        # Upload
        #----------------------------------------------------------------------------
        if _Changed and upload:
            djOrg.save()
            
        
    #------------------------------------------------------------------  
    def handle(self, *args, **options):
    #------------------------------------------------------------------  

        if options['debug']:
            logger.setLevel(logging.DEBUG)
            verbose=2
        else:
            logger.setLevel(logging.INFO)
            verbose=0
            
        # qryOrg = Organism.objects.all()

        if options['orgid']:
            logger.info(f" [{self.command_process}] Processing OrgID {options['orgid']} ...")
            self._process_organism(options['orgid'], upload=options['upload'],overwrite=options['overwrite'],verbose=verbose)

        elif options['all']:
            qryOrg = Organism.objects.all().values_list('organism_id')
            nOrg = qryOrg.count()
            logger.info(f" [{self.command_process}] Processing All {nOrg} Organisms  ...")

            for _org_id in tqdm(qryOrg, total=nOrg, desc='Organism'):
                self._process_organism(_org_id[0], upload=options['upload'],overwrite=options['overwrite'],verbose=verbose)


        
        
                    
