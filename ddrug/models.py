from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import *

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError

from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary
from dorganism.models import Organism, Organism_Batch
from dscreen.models import Screen_Run
#-------------------------------------------------------------------------------------------------
# Drugs related Application Model
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Drug(AuditModel):
    """
    List of Drugs, DrugCombinations, DrugScreens 
    """
#=================================================================================================
    HEADER_FIELDS = {
        "drug_id":{"Drug ID": {"drug_id": LinkList["drug_id"] } },
        "drug_name":"Drug Name",
        "drug_othernames":"Other Names",
        "drug_codes":"Drug Codes",
        "drug_type":"Drug Type",
        "drug_class":"Drug Class",
        "drug_subclass":"Sub Class",
    }
    
    Choice_Dictionary = {
        'drug_type':'Drug_Type',
        'max_phase':'Max_Phase',
    }

    CARDS_FIELDS= {
        "antimicro" : "MIC",

    }

  
    ID_SEQUENCE = 'Drug'
    ID_PREFIX = 'AMD'
    ID_PAD = 5

    drug_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Drug ID")
    drug_name = models.CharField(max_length=50, unique=True, verbose_name = "Drug Name")
    drug_othernames = ArrayField(models.CharField(max_length=60, blank=True),size=30, null=True,verbose_name = "Other Names")
    drug_codes = ArrayField(models.CharField(max_length=10, blank=True),size=30, null=True)
    drug_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Drug Type", on_delete=models.DO_NOTHING,
        db_column="drug_type", related_name="%(class)s_drugtype")
    drug_note = models.CharField(max_length=50, blank=True, verbose_name = "Drug Note")
    drug_panel=ArrayField(models.CharField(max_length=20, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)

    drug_target = models.CharField(max_length=50, blank=True,  verbose_name = "Target")
    drug_subtarget = models.CharField (max_length=50, blank=True,  verbose_name = "SubTarget")
    moa = models.CharField(blank=True, max_length=50, verbose_name = "MoA")
    drug_class = models.CharField(max_length=50, blank=True,  verbose_name = "Class")
    drug_subclass = models.CharField(max_length=100, blank=True, verbose_name = "SubClass")

    antimicro = models.CharField(max_length=25, blank=True, verbose_name = "Antimicro")
    antimicro_class = models.CharField(max_length=80, blank=True, verbose_name = "Antimicro Type")

    max_phase = models.CharField(max_length=2, default="0", blank=True, verbose_name = "max Approval")
    approval_note = models.CharField(max_length=50, blank=True, verbose_name = "Approval Note")
    admin_routes = models.CharField(max_length=50, blank=True, verbose_name = "Administrations")
    application =  models.CharField(max_length=50, blank=True, verbose_name = "Application")

    n_compounds = models.IntegerField(default=0, blank=True, verbose_name = "#Cmpds")	
    chembl =    models.CharField(max_length=15, blank=True, verbose_name = "ChEMBL")	
    drugbank =  models.CharField(max_length=10, blank=True, verbose_name = "DrugBank")	
    cas	=       models.CharField(max_length=15, blank=True, verbose_name = "CAS")	
    pubchem	=   models.CharField(max_length=15, blank=True, verbose_name = "PubChem")	
    chemspider= models.CharField(max_length=15, blank=True, verbose_name = "ChemSpider")	
    unii = 	    models.CharField(max_length=12, blank=True, verbose_name = "UNII")	
    kegg = 	    models.CharField(max_length=10, blank=True, verbose_name = "KEGG")	
    comptox	=   models.CharField(max_length=20, blank=True, verbose_name = "CompTox")	
    echa = 	    models.CharField(max_length=15, blank=True, verbose_name = "ECHA")	
    chebi =     models.CharField(max_length=15, blank=True, verbose_name = "ChEBI")	
    uq_imb =    models.CharField(max_length=15, blank=True, verbose_name = "IMB")	
    vendor =    models.CharField(max_length=15, blank=True, verbose_name = "Vendor")	
    vendor_catno = models.CharField(max_length=15, blank=True, verbose_name = "CatNo")	
    mw = models.DecimalField(default=0, max_digits=12, decimal_places=2, blank=True, verbose_name = "MW")	
    mf = models.CharField(max_length=25, blank=True, verbose_name = "MF")	
    smiles = models.CharField(max_length=2048, blank=True, verbose_name = "SMILES")

    # Require RDKit-PostgreSQL 
    smol = models.MolField(blank=True, null=True, verbose_name = "MOL")	
    torsionbv = models.BfpField(null=True)	
    ffp2 = models.BfpField(null=True, verbose_name = "FFP2")
    mfp2 = models.BfpField(null=True, verbose_name = "MFP2")
    salt_form = models.CharField(blank=True, max_length=15, verbose_name = "SaltForm")	

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'drug'
        ordering=['drug_name']
        indexes = [
            models.Index(name="drug_dname_idx", fields=['drug_name']),
            #GistIndex(name="drug_smol_idx",fields=['smol']),
            #GistIndex(name="drug_ffp2_idx",fields=['ffp2']),
            #GistIndex(name="drug_mfp2_idx",fields=['mfp2'])
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.drug_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.drug_name} ({self.drug_id})"

#    #------------------------------------------------
#     @classmethod
#     def str_DrugID(cls,DrugNo) -> str:
#         return(f"{cls.ID_PREFIX}{DrugNo:05d}")

#     #------------------------------------------------
#     @classmethod
#     def find_Next_DrugID(cls) -> str:
#         Drug_IDSq=Sequence(cls.ID_SEQUENCE)
#         Drug_nextID = next(Drug_IDSq)
#         Drug_strID = cls.str_DrugID(Drug_nextID)
#         while cls.exists(None,Drug_strID):
#             Drug_nextID = next(Drug_IDSq)
#             Drug_strID = cls.str_DrugID(Drug_nextID)
#         return(Drug_strID)    

    #------------------------------------------------
    @classmethod
    def get(cls,DrugName,DrugID=None,verbose=0):
    # Returns an instance by drug_name or durg_id
        try:
            if DrugName:
                retInstance = cls.objects.get(drug_name=DrugName)
            elif DrugID:
                retInstance = cls.objects.get(drug_id=DrugID)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if DrugName:
                    print(f"[Drug Not Found] {DrugName} ")
                elif DrugID:
                    print(f"[Drug Not Found] {DrugID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,DrugName,DrugID=None,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        if DrugName:
            retValue = cls.objects.filter(drug_name=DrugName).exists()
        elif DrugID:
            retValue = cls.objects.filter(drug_id=DrugID).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    @classmethod
    def update_all_fp(cls):
        pass
        # Require RDKit-PostgreSQL
        # cls.objects.update(ffp2=FEATMORGANBV_FP('smol'),
        #                    mfp2=MORGANBV_FP('smol'), 
        #                    torsionbv=TORSIONBV_FP('smol')
        #                    )

    #------------------------------------------------
    @classmethod
    def smiles2mol(cls,Smiles,verbose=0):
        try:
            xmol = Chem.MolFromSiles(Smiles)
        except:
            xmol = None
            if verbose:
                print(f"[Invalid SMILES] {Smiles} ")
        return(xmol)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.drug_id:
            self.drug_id = self.next_id()
            if self.drug_id: 
                super(Drug, self).save(*args, **kwargs)
                # self.__dict__.update(ffp2=FEATMORGANBV_FP('smol'), mfp2=MORGANBV_FP('smol'), torsionbv=TORSIONBV_FP('smol'))

                # Require RDKit-PostgreSQL
                # Drug.objects.filter(drug_id=self.drug_id).update(
                #     ffp2=FEATMORGANBV_FP('smol'), 
                #     mfp2=MORGANBV_FP('smol'), 
                #     torsionbv=TORSIONBV_FP('smol')
                #     )
        else:
            # Require RDKit-PostgreSQL
            # self.__dict__.update(
            #     ffp2=FEATMORGANBV_FP('smol'), 
            #     mfp2=MORGANBV_FP('smol'), 
            #     torsionbv=TORSIONBV_FP('smol')
            #     )
            super(Drug, self).save(*args, **kwargs) 

#=================================================================================================
class Breakpoint(AuditModel):
    """
    List of Breakpoints 
    """
#=================================================================================================
    HEADER_FIELDS = {
       'drug_id.drug_name':{'Drug Name': {'drug_id.drug_id':LinkList["drug_id"]}},
       'org_name':'org_name', 
       'org_rank':'org_rank', 
       'notorg_name':'notorg_name', 
       'notorg_rank':'notorg_rank',
       'med_application':'med_application',
        'bp_type':'bp_type', 
        'bp_res_gt':'bp_res_gt', 
        'bp_sens_le':'bp_sens_le',
        'bp_unit':'bp_unit', 
        'bp_comb':'bp_comb', 
        'bp_source':'bp_source', 
        'bp_source_version':'bp_source_version',

    }

    Choice_Dictionary= {
        'org_rank':'Tax_Rank',
        'notorg_rank':'Tax_Rank',
        'bp_type':'BP_Type',
    }

    drug_id = models.ForeignKey(Drug, null=False, blank=False, verbose_name = "Drug ID", on_delete=models.DO_NOTHING,
        db_column="drug_id", related_name="%(class)s_drug_id")
    org_name = models.CharField(max_length=50, blank=True, verbose_name = "Organism") 
    org_rank = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Rank", on_delete=models.DO_NOTHING,
        db_column="org_rank", related_name="%(class)s_orgtype")
    notorg_name = models.CharField(max_length=50, blank=True, verbose_name = "Not(Organism)") 
    notorg_rank = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Not(Rank)", on_delete=models.DO_NOTHING,
        db_column="notorg_rank", related_name="%(class)s_notorgtype")
    med_application = models.CharField(max_length=50, blank=True, verbose_name = "Application")
    bp_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "BP Type", on_delete=models.DO_NOTHING,
        db_column="bp_type", related_name="%(class)s_bptype")
    bp_res_gt = models.DecimalField(max_digits=9, decimal_places=3, blank=False, verbose_name = ">Res") 
    bp_sens_le = models.DecimalField(max_digits=9, decimal_places=3, blank=False, verbose_name = "<=Sens") 
    bp_unit = models.CharField(max_length=5, blank=False, verbose_name = "Unit") 
    bp_comb = models.CharField(max_length=20, blank=True, verbose_name = "Combination") 
    bp_source = models.CharField(max_length=50, blank=False, verbose_name = "BP Source") 
    bp_source_version = models.CharField(max_length=50, blank=False, verbose_name = "BP Version") 

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'breakpoint'
        ordering=['drug_id','org_rank','org_name','bp_type']
        indexes = [
            models.Index(name="bp_drug_idx", fields=['drug_id']),
            models.Index(name="bp_org_idx",fields=['org_name']),
            models.Index(name="bp_nrnk_idx",fields=['notorg_rank']),
            models.Index(name="bp_norg_idx",fields=['notorg_name']),
            models.Index(name="bp_ornk_idx",fields=['org_rank']),
            models.Index(name="bp_src_idx",fields=['bp_source']),
            models.Index(name="bp_btyp_idx",fields=['bp_type']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.drug_id} : {self.org_name} ({self.org_rank}) Not: {self.notorg_name} ({self.notorg_rank}) Med: {self.med_application} BP: {self.bp_type} {self.bp_source}"

    #------------------------------------------------
    def info(self) -> str:
        if self.org_rank:
            return(f"{self.bp_source} ({self.org_rank})")
        elif self.notorg_rank:
            return(f"{self.bp_source} (Not({self.notorg_rank}))")
        else:
            return(f"{self.bp_source}")

   #------------------------------------------------
    @classmethod
    def get(cls,DrugID, OrgName, OrgRank, NotOrgName, NotOrgRank, MedAppl, BPType, BPSource, verbose=0):
    # Returns an instance if found 
        try:
            retInstance = cls.objects.get(drug_id=DrugID, 
                                          org_name=OrgName, org_rank=OrgRank, notorg_name=NotOrgName, notorg_rank=NotOrgRank, med_application = MedAppl,  
                                          bp_type=BPType, bp_source=BPSource)
        except:
            if verbose:
                print(f"[Breakpoint  Not Found] {DrugID} {OrgName} {OrgRank} Not: {NotOrgName} {NotOrgRank} Med: {MedAppl} BP: {BPType} {BPSource}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,DrugID, OrgName, OrgRank, NotOrgName, NotOrgRank, MedAppl, BPType, BPSource, verbose=0):
    # Returns if an instance exists 
        return cls.objects.filter(drug_id=DrugID, 
                                          org_name=OrgName, org_rank=OrgRank, notorg_name=NotOrgName, notorg_rank=NotOrgRank, med_application = MedAppl,  
                                          bp_type=BPType, bp_source=BPSource).exists()


#=================================================================================================
class VITEK_Card(AuditModel):
#     """
#     List of VITEK Cards
#     """
#=================================================================================================
    HEADER_FIELDS = {
        "orgbatch_id.organism_id.organism_id":{'Organism ID': {'orgbatch_id.organism_id.organism_id':LinkList["organism_id"]}},
        "orgbatch_id.batch_id":"Batch",
        #"orgbatch_id":"Orgbatch",
        "card_barcode":"Barcode",
        "card_type":"Card Type",
        "card_code":"Card Code",
        "expiry_date":"Expiry",
        "proc_date":"Processed",
        "analysis_time":"Time",
        "instrument":"Instrument",
    }

    Choice_Dictionary= {
        'card_type':'Card_Type',
    }

    card_barcode = models.CharField(max_length=25, primary_key=True, verbose_name = "Card Barcode") 
    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    card_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Card Type", on_delete=models.DO_NOTHING,
        db_column="card_type", related_name="%(class)s_cardtype")
    card_code = models.CharField(max_length=15, blank=False, verbose_name = "Card Code") 
    expiry_date = models.DateField(blank=True, verbose_name = "Expiry Date")
    instrument = models.CharField(max_length=50, blank=True, verbose_name = "Instrument") 
    proc_date = models.DateField(blank=False, verbose_name = "Processing Date")
    analysis_time = models.CharField(max_length=15, blank=True, verbose_name = "Analysis Time") 

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'vitek_card'
        ordering=['card_type','card_code','card_barcode']
        indexes = [
            models.Index(name="vcard_ctype_idx",fields=['card_type']),
            models.Index(name="vcard_ccode_idx",fields=['card_code']),
            models.Index(name="vcard_bcode_idx",fields=['card_barcode']),
            models.Index(name="vcard_orgbatch_idx",fields=['orgbatch_id']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.card_barcode}"
    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.orgbatch_id} {self.card_code} {self.card_barcode}"

   #------------------------------------------------
    @classmethod
    def get(cls,CardBarcode,verbose=0):
    # Returns an instance if found by Card Barcode
        try:
            retInstance = cls.objects.get(card_barcode=CardBarcode)
        except:
            if verbose:
                print(f"[Vitek Card Not Found] {CardBarcode}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,CardBarcode,verbose=0):
    # Returns if an instance exists by Card Barcode
        return cls.objects.filter(card_barcode=CardBarcode).exists()


#=================================================================================================
class VITEK_AST(AuditModel):
    """
      Antimicrobial Suceptibility Testing (AST) data from VITEK Cards
    """
#=================================================================================================
    HEADER_FIELDS = {
        "card_barcode.orgbatch_id.organism_id.organism_id":{'Organism ID': {'card_barcode.orgbatch_id.organism_id.organism_id':LinkList["organism_id"]}},
        "card_barcode.orgbatch_id.batch_id":"Batch",
        "card_barcode.orgbatch_id.organism_id.organism_name":"Organism Name",
        "drug_id.drug_name":{'Drug Name': {'drug_id.drug_id':LinkList['drug_id']}},
        "drug_id.drug_codes":"Codes",
        "mic":"MIC",
        "bp_profile":"BP",
        "organism":"Selected Organism",
        "bp_comment":"Comment",
        "bp_source":"Source",
    }

    card_barcode = models.ForeignKey(VITEK_Card, null=False, blank=False, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING,
        db_column="card_barcode", related_name="%(class)s_card_barcode") 
    drug_id = models.ForeignKey(Drug, null=False, blank=False, verbose_name = "Drug ID", on_delete=models.DO_NOTHING,
        db_column="drug_id", related_name="%(class)s_drug_id")
    mic = models.CharField(max_length=50, blank=True, verbose_name = "MIC")
    process = models.CharField(max_length=50, blank=True, verbose_name = "Vitek Process")
    bp_profile = models.CharField(max_length=5, blank=True, verbose_name = "Break Point")
    bp_comment = models.CharField(max_length=120, blank=True, verbose_name = "Comment")
    bp_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")
    selection = models.CharField(max_length=20, blank=True, verbose_name = "Selection")
    organism = models.CharField(max_length=120, blank=True, verbose_name = "Selected Organism")
    filename = models.CharField(max_length=120, blank=True, verbose_name = "PDF Filename")
    page_no = models.IntegerField(default=0, blank=True, verbose_name = "PDF PageNo")

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'vitek_ast'
        ordering=['drug_id','bp_source','card_barcode']
        indexes = [
            models.Index(name="vast_drugid_idx",fields=['drug_id']),
            models.Index(name="vast_mic_idx",fields=['mic']),
            models.Index(name="vast_bprofile_idx",fields=['bp_profile']),
            models.Index(name="vast_bpsource_idx",fields=['bp_source']),
            models.Index(name="vast_fname_idx",fields=['filename']),
            models.Index(name="vast_orgname_idx",fields=['organism']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = ""
        if self.drug_id:    
            if self.drug_id is not None:
                retStr += f"{self.drug_id.drug_name} "
            else:
                retStr += f"{self.drug_id} "
        if self.card_barcode:
            if self.card_barcode is not None:
                retStr += f"{self.card_barcode.orgbatch_id} {self.card_barcode.card_code}"
            else:
                retStr += f"{self.card_barcode} "
        return(retStr)

    #------------------------------------------------
    @classmethod
    def get(cls,CardBarcode,DrugID,Source,OrgName,verbose=0):
    # Returns an instance if found by (CardBarcode,DrugID,Source)
        try:
            retInstance = cls.objects.get(card_barcode=CardBarcode,drug_id=DrugID,bp_source=Source,organism=OrgName)
        except:
            if verbose:
                print(f"[Vitek AST Not Found] {CardBarcode} {DrugID} {Source} {OrgName}")
            retInstance = None
        return(retInstance)
    #------------------------------------------------
    @classmethod
    def exists(cls,CardBarcode,DrugID,Source,OrgName,verbose=0):
    # Returns an instance if found by (CardBarcode,DrugID,Source)
        return cls.objects.filter(card_barcode=CardBarcode,drug_id=DrugID,bp_source=Source,organism=OrgName).exists()
    
   


    
#=================================================================================================
class VITEK_ID(AuditModel):
    """
      Identification Testing (ID) data from VITEK Cards
    """
#=================================================================================================
    HEADER_FIELDS = {
        "card_barcode.orgbatch_id.organism_id.organism_id":{'Org ID': {'card_barcode.orgbatch_id.organism_id.organism_id':LinkList["organism_id"]}},
        "card_barcode.orgbatch_id.batch_id":"Batch",
        "card_barcode.orgbatch_id.organism_id.organism_name":"Organism",
        "id_organism":"Identification",
        "id_probability":"Probability",
        "id_confidence":"Confidence",
        "id_source":"Source",
        "process":"Vitek Process",
    }

    card_barcode = models.ForeignKey(VITEK_Card, null=False, blank=False, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING,
        db_column="card_barcode", related_name="%(class)s_card_barcode") 
    process = models.CharField(max_length=50, blank=True, verbose_name = "Vitek Process")
    id_organism = models.CharField(max_length=120,  blank=True, verbose_name = "ID Organism")
    id_probability = models.CharField(max_length=120, blank=True,  verbose_name = "ID Probability")
    id_confidence = models.CharField(max_length=120, blank=True,  verbose_name = "ID Confidence")
    id_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")
    filename = models.CharField(max_length=120, blank=True, verbose_name = "PDF Filename")
    page_no = models.IntegerField(default=0, blank=True, verbose_name = "PDF PageNo")  

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'vitek_id'
        ordering=['card_barcode']
        indexes = [
            models.Index(name="vid_barcode_idx",fields=['card_barcode']),
            models.Index(name="vid_idorg_idx",fields=['id_organism']),
            models.Index(name="vid_idprob_idx",fields=['id_probability']),
            models.Index(name="vid_idconf_idx",fields=['id_confidence']),
            models.Index(name="vid_fname_idx",fields=['filename']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = "<empty>"
        if self.card_barcode:
            if self.card_barcode is not None:
                retStr = f"{self.id_organism} {self.card_barcode.orgbatch_id} {self.card_barcode.card_code}"
            else:
                retStr = f"{self.card_barcode}"
        return retStr

   #------------------------------------------------
    @classmethod
    def get(cls,CardBarcode,verbose=0):
    # Returns an instance if found by CardBarcode
        try:
            retInstance = cls.objects.get(card_barcode=CardBarcode)
        except:
            if verbose:
                print(f"[Vitek AST Not Found] {CardBarcode} ")
            retInstance = None
        return(retInstance)
   #------------------------------------------------
    @classmethod
    def exists(cls,CardBarcode,verbose=0):
    # Returns if an instance exists by CardBarcode
        return cls.objects.filter(card_barcode=CardBarcode).exists()


#=================================================================================================
class MIC_COADD(AuditModel):
    """
     Antibiogram from CO-ADD screening    
    """
#=================================================================================================
    HEADER_FIELDS = {
        "orgbatch_id.organism_id.organism_id":{'Organism ID': {'orgbatch_id.organism_id.organism_id':LinkList["organism_id"]}},
        "orgbatch_id.batch_id":"Batch",
        "orgbatch_id.organism_id.organism_name":"Organism Name",
        "drug_id.drug_name":{'Drug Name': {'drug_id.drug_id':LinkList['drug_id']}},
        "bp_profile":"Break Point",
        #"mic_type":"Type",
        "mic":"MIC",
        "run_id":"Run ID",
        #"media":"Media",
    }
    
    Choice_Dictionary = {
        'mic_type':'MIC_Type',
        'plate_size':'Plate_Size',
        'plate_material':'Plate_Material',
        'media':'Media',
    }

    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id") 
    drug_id = models.ForeignKey(Drug, null=False, blank=False, verbose_name = "Drug ID", on_delete=models.DO_NOTHING,
        db_column="drug_id", related_name="%(class)s_drug_id")
    
    mic = models.CharField(max_length=50, blank=True, verbose_name = "MIC")
    mic_unit = models.CharField(max_length=20, blank=True, verbose_name = "Unit")
    mic_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MIC Type", on_delete=models.DO_NOTHING,
         db_column="mic_type", related_name="%(class)s_mictype")
    bp_profile = models.CharField(max_length=5, blank=True, verbose_name = "Break Point")
    bp_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")

    run_id = models.ForeignKey(Screen_Run, null=False, blank=False, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
        db_column="run_id", related_name="%(class)s_run_id") 
    # run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID")
    testplate_id = models.CharField(max_length=25, blank=True, verbose_name = "PlateID")
    testwell_id = models.CharField(max_length=5, blank=True, verbose_name = "WellID")

    plate_size = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Size", on_delete=models.DO_NOTHING,
        db_column="plate_size", related_name="%(class)s_platesize")
    plate_material = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Plate Material", on_delete=models.DO_NOTHING,
        db_column="plate_material", related_name="%(class)s_material")

    # Possible update to ForeignKey (JZG) 
    #media = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Media", on_delete=models.DO_NOTHING,
    #    db_column="media", related_name="%(class)s_Media+")
    media = models.CharField(max_length=40, blank=True, verbose_name = "Media")
    dye = models.CharField(max_length=40, blank=True, verbose_name = "Dye")

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'mic_coadd'
        ordering=['drug_id','orgbatch_id']
        indexes = [
             models.Index(name="micc_drugid_idx",fields=['drug_id']),
             models.Index(name="micc_mic_idx",fields=['mic']),
             models.Index(name="micc_bprofile_idx",fields=['bp_profile']),
             models.Index(name="micc_bpsource_idx",fields=['bp_source']),
            #  models.Index(name="micc_source_idx",fields=['source']),
             models.Index(name="micc_size_idx",fields=['plate_size']),
             models.Index(name="micc_material_idx",fields=['plate_material']),
             models.Index(name="micc_media_idx",fields=['media']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = f"{self.drug_id} {self.orgbatch_id}"
        return(retStr)

    def __repr__(self) -> str:
        retStr = f"{self.drug_id.drug_name} {self.orgbatch_id} {self.mic} {str(self.run_id)}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,DrugID,TestPlateID,TestWellID,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID,drug_id=DrugID,testplate_id=TestPlateID,testwell_id=TestWellID)
        except:
            if verbose:
                print(f"[MIC Not Found] {OrgBatchID} {DrugID} {TestPlateID} {TestWellID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,DrugID,TestPlateID,TestWellID,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        return cls.objects.filter(orgbatch_id=OrgBatchID,drug_id=DrugID,testplate_id=TestPlateID,testwell_id=TestWellID).exists()


#=================================================================================================
class MIC_Pub(AuditModel):
    """
     Antibiogram from Public sources    
    """
#=================================================================================================
    HEADER_FIELDS   = {   
        "organism_id.organism_id":{'Organism ID': {'organism_id.organism_id':LinkList['organism_id']}}, 
        "organism_id.organism_name":"Organism",
        "drug_id.drug_name":{'Drug Name': {'drug_id.drug_id':LinkList['drug_id']}},
        "bp_profile":"BP",
        "mic":"MIC",
        "zone_diameter": "Zone",
        "mic_type":"Type",
        "source":"Source",
    }
    Choice_Dictionary = {
        'mic_type':'MIC_Type',
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="organism_id", related_name="%(class)s_organism_id") 
    drug_id = models.ForeignKey(Drug, null=False, blank=False, verbose_name = "Drug ID", on_delete=models.DO_NOTHING,
        db_column="drug_id", related_name="%(class)s_drug_id")
    
    mic = models.CharField(max_length=50, blank=True, verbose_name = "MIC")
    mic_unit = models.CharField(max_length=20, blank=True, verbose_name = "Unit")
    zone_diameter = models.CharField(max_length=20, blank=True, verbose_name = "Zone Diameter")
    mic_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "MIC Type", on_delete=models.DO_NOTHING,
         db_column="mic_type", related_name="%(class)s_mictype")
    source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
    bp_profile = models.CharField(max_length=5, blank=True, verbose_name = "Break Point")
    bp_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")


    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'mic_pub'
        ordering=['drug_id','organism_id']
        indexes = [
             models.Index(name="micp_drugid_idx",fields=['drug_id']),
             models.Index(name="micp_mic_idx",fields=['mic']),
             models.Index(name="micp_bprofile_idx",fields=['bp_profile']),
             models.Index(name="micp_bpsource_idx",fields=['bp_source']),
             models.Index(name="micp_source_idx",fields=['source']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        retStr = ""
        if self.drug_id:    
            if self.drug_id is not None:
                retStr += f"{self.drug_id.drug_name} "
            else:
                retStr += f"{self.drug_id} "
        retStr += f"{self.organism_id} {self.mic} {self.source}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,OrgID,DrugID,Source,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        try:
            retInstance = cls.objects.get(organism_id=OrgID,drug_id=DrugID,source=Source)
        except:
            if verbose:
                print(f"[MIC Not Found] {OrgID} {DrugID} {Source}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgID,DrugID,Source,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        return cls.objects.filter(organism_id=OrgID,drug_id=DrugID,source=Source).exists()

   
