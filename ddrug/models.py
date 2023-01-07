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

#-------------------------------------------------------------------------------------------------
# Organism Drugs Application Model
#-------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------
class Drug(AuditModel):
    """
    List of Drugs, DrugCombinations, DrugScreens 
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary = {
        'drug_type':'Drug_Type',
        'max_phase':'Max_Phase',
    }

    drug_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Drug ID")
    drug_name = models.CharField(max_length=50, unique=True, verbose_name = "Drug Name")
    drug_othernames = ArrayField(models.CharField(max_length=60, blank=True),size=30, null=True)
    drug_codes = ArrayField(models.CharField(max_length=10, blank=True),size=30, null=True)
    drug_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Drug Type", on_delete=models.DO_NOTHING,
        db_column="drug_type", related_name="%(class)s_DrugType+")
    drug_note = models.CharField(max_length=50, blank=True, verbose_name = "Drug Notes")
    drug_panel=ArrayField(models.CharField(max_length=20, null=True, blank=True), size=20, verbose_name = "Panel", null=True, blank=True)

    drug_target = models.CharField(max_length=50, blank=True,  verbose_name = "Target")
    drug_subtarget = models.CharField (max_length=50, blank=True,  verbose_name = "SubTarget")
    moa = models.CharField(blank=True, max_length=50, verbose_name = "MoA")
    drug_class = models.CharField(max_length=50, blank=True,  verbose_name = "Class")
    drug_subclass = models.CharField(max_length=50, blank=True, verbose_name = "SubClass")

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
    smol = models.MolField(blank=True, null=True, verbose_name = "MOL")	
    ffp2 = models.BfpField(null=True, verbose_name = "FFP2")
    #salt_form = models.CharField(blank=True, max_length=15, verbose_name = "SaltForm")	

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'drug'
        ordering=['drug_name']
        indexes = [
            models.Index(name="drug_dname_idx", fields=['drug_name']),
            GistIndex(name="drug_smol_idx",fields=['smol']),
            GistIndex(name="drug_ffp2_idx",fields=['ffp2'])
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.drug_name} ({self.drug_id})"

    #------------------------------------------------
    @classmethod
    def exists(self,DrugName,DrugID=None,verbose=0):
    #
    # Returns an instance by drug_name or durg_id
    #
        try:
            if DrugName:
                retInstance = self.objects.get(drug_name=DrugName)
            elif DrugID:
                retInstance = self.objects.get(drug_id=DrugID)
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
    def update_all_fp(self):
        self.objects.update(ffp2=FEATMORGANBV_FP('smol'))

    #------------------------------------------------
    @classmethod
    def smiles2mol(self,Smiles,verbose=0):
        try:
            xmol = Chem.MolFromSiles(Smiles)
        except:
            xmol = None
            if verbose:
                print(f"[Invalid SMILES] {Smiles} ")
        return(xmol)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        self.__dict__.update(ffp2=FEATMORGANBV_FP('smol'))
        super(Drug, self).save(*args, **kwargs)
            
    # -------------------------------------------------
    def get_values(self, fields=DRUG_FIELDs):
        value_list=super(Drug, self).get_values(fields)
        return value_list
        

#-------------------------------------------------------------------------------------------------
class VITEK_Card(AuditModel):
#     """
#     List of VITEK Cards
#     """
# #-------------------------------------------------------------------------------------------------
    Choice_Dictionary= {
        'card_type':'Card_Type',
    }

    card_barcode = models.CharField(max_length=25, primary_key=True, verbose_name = "Card Barcode") 
    orgbatch_id = models.ForeignKey(Organism_Batch, null=False, blank=False, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING,
        db_column="orgbatch_id", related_name="%(class)s_orgbatch_id+") 
    card_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Card Type", on_delete=models.DO_NOTHING,
        db_column="card_type", related_name="%(class)s_CardType+")
    card_code = models.CharField(max_length=15, blank=False, verbose_name = "Card Code") 
    expiry_date = models.DateField(blank=True, verbose_name = "Expiry Date")
    instrument = models.CharField(max_length=50, blank=True, verbose_name = "Vitek SN") 
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
        return f"{self.orgbatch_id} {self.card_type} {self.card_code}"

   #------------------------------------------------
    @classmethod
    def exists(self,CardBarcode,verbose=0):
    #
    # Returns an instance if found by Card Barcode
    #
        try:
            retInstance = self.objects.get(card_barcode=CardBarcode)
        except:
            if verbose:
                print(f"[Vitek Card Not Found] {CardBarcode}")
            retInstance = None
        return(retInstance)


#-------------------------------------------------------------------------------------------------
class VITEK_AST(AuditModel):
#     """
#      Antimicrobial Susceptibility Testing (AST) data from VITEK Cards
#     """
# #-------------------------------------------------------------------------------------------------

    card_barcode = models.ForeignKey(VITEK_Card, null=False, blank=False, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING,
        db_column="card_barcode", related_name="%(class)s_card_barcode+") 
    drug_id = models.ForeignKey(Drug, null=False, blank=False, verbose_name = "Drug ID", on_delete=models.DO_NOTHING,
        db_column="drug_id", related_name="%(class)s_drug_id+")
    mic = models.CharField(max_length=50, blank=True, verbose_name = "MIC")
    process = models.CharField(max_length=50, blank=True, verbose_name = "Vitek Process")
    bp_profile = models.CharField(max_length=5, blank=True, verbose_name = "Break Point")
    bp_comment = models.CharField(max_length=120, blank=True, verbose_name = "Comment")
    bp_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")
    selection = models.CharField(max_length=20, blank=True, verbose_name = "Selection")
    organism = models.CharField(max_length=120, blank=True, verbose_name = "Organism")
    filename = models.CharField(max_length=120, blank=True, verbose_name = "PDF Filename")
    page_no = models.IntegerField(default=0, blank=True, verbose_name = "PDF PageNo")  

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'vitek_ast'
        ordering=['drug_id','card_barcode']
        indexes = [
            models.Index(name="vast_drugid_idx",fields=['drug_id']),
            models.Index(name="vast_mic_idx",fields=['mic']),
            models.Index(name="vast_bprofile_idx",fields=['bp_profile']),
            models.Index(name="vast_fname_idx",fields=['filename']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.drug_id} {self.card_barcode.orgbatch_id} {self.card_barcode.card_code}"

   #------------------------------------------------
    @classmethod
    def exists(self,CardBarcode,DrugID,verbose=0):
    #
    # Returns an instance if found by CardBarcode and DrugID
    #
        try:
            retInstance = self.objects.get(card_barcode=CardBarcode,drug_id=DrugID)
        except:
            if verbose:
                print(f"[Vitek AST Not Found] {CardBarcode} {DrugID}")
            retInstance = None
        return(retInstance)

# #-------------------------------------------------------------------------------------------------
class VITEK_ID(AuditModel):
#     """
#      Identification Testing (ID) data from VITEK Cards
#     """
# #-------------------------------------------------------------------------------------------------

    card_barcode = models.ForeignKey(VITEK_Card, null=False, blank=False, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING,
        db_column="card_barcode", related_name="%(class)s_card_barcode+") 
    process = models.CharField(max_length=50, blank=True, verbose_name = "Vitek Process")
    id_organism = models.CharField(max_length=120,  blank=True, verbose_name = "ID Organism")
    id_probability = models.CharField(max_length=120, blank=True,  verbose_name = "ID Organism")
    id_confidence = models.CharField(max_length=120, blank=True,  verbose_name = "ID Organism")
    id_source = models.CharField(max_length=20,  blank=True, verbose_name = "Source")
    filename = models.CharField(max_length=120, blank=True, verbose_name = "PDF Filename")
    page_no = models.IntegerField(default=0, blank=True, verbose_name = "PDF PageNo")  

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'vitek_id'
        #ordering=['card_barcode']
        indexes = [
            models.Index(name="vid_idorg_idx",fields=['id_organism']),
            models.Index(name="vid_idprob_idx",fields=['id_probability']),
            models.Index(name="vid_idconf_idx",fields=['id_confidence']),
            models.Index(name="vid_fname_idx",fields=['filename']),
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.id_organism} {self.card_barcode.orgbatch_id} {self.card_barcode.card_code}"

   #------------------------------------------------
    @classmethod
    def exists(self,CardBarcode,verbose=0):
    #
    # Returns an instance if found by CardBarcode and DrugID
    #
        try:
            retInstance = self.objects.get(card_barcode=CardBarcode)
        except:
            if verbose:
                print(f"[Vitek AST Not Found] {CardBarcode} ")
            retInstance = None
        return(retInstance)

#-------------------------------------------------------------------------------------------------
