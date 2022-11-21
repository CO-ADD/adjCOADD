from model_utils import Choices
from sequences import Sequence
from django_rdkit import models
from django.contrib.postgres.fields import ArrayField
from django.db import transaction, IntegrityError

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
        'Drug_Type':'Drug_Type',
        'Max_Phase':'Max_Phase',
    }

    Drug_ID = models.CharField(primary_key=True, max_length=15, verbose_name = "Drug ID")
    Drug_Name = models.CharField(blank=False, null=False, unique=True, max_length=50, verbose_name = "Drug Name")
    Drug_OtherName = models.CharField(blank=True, null=True, max_length=150, verbose_name = "Drug OtherName")
    Drug_Code = models.CharField(blank=True, null=True, unique=True, max_length=10, verbose_name = "Drug Code")
    Drug_Type = models.CharField(blank=True, null=True, max_length=15, verbose_name = "Drug Type")
    Drug_Note = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Drug Notes")

    ATC_CODE = models.CharField(blank=True, null=True, max_length=50, verbose_name = "ATC Code")
    Drug_Target = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Target")
    Drug_SubTarget = models.CharField(blank=True, null=True, max_length=50, verbose_name = "SubTraget")
    Drug_Class = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Class")
    Drug_SubClass = models.CharField(blank=True, null=True, max_length=50, verbose_name = "SubClass")

    Mode_Action = models.CharField(blank=True, null=True, max_length=50, verbose_name = "MoA")
    Antimicro = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Antimicro")
    Antimicro_Type = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Antimicro Type")

    Max_Phase = models.CharField(blank=True, null=True, max_length=2, verbose_name = "max Approval")
    Approval_Note = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Approval Notes")
    Admin_Route = models.CharField(blank=True, null=True, max_length=50, verbose_name = "Administration")
    Application =  models.CharField(blank=True, null=True, max_length=50, verbose_name = "Application")

    N_Compounds = models.IntegerField(blank=True, null=True, verbose_name = "#Cmpds")	
    ChEMBL =    models.CharField(blank=True, null=True, max_length=15, verbose_name = "ChEMBL")	
    DrugBank =  models.CharField(blank=True, null=True, max_length=15, verbose_name = "DrugBank")	
    CAS	=       models.CharField(blank=True, null=True, max_length=15, verbose_name = "CAS")	
    PubChem	=   models.CharField(blank=True, null=True, max_length=15, verbose_name = "PubChem")	
    ChemSpider= models.CharField(blank=True, null=True, max_length=15, verbose_name = "ChemSpider")	
    UNII = 	    models.CharField(blank=True, null=True, max_length=15, verbose_name = "UNII")	
    KEGG = 	    models.CharField(blank=True, null=True, max_length=15, verbose_name = "KEGG")	
    CompTox	=   models.CharField(blank=True, null=True, max_length=15, verbose_name = "CompTox")	
    ECHA = 	    models.CharField(blank=True, null=True, max_length=15, verbose_name = "ECHA")	
    ChEBI =     models.CharField(blank=True, null=True, max_length=15, verbose_name = "ChEBI")	
    UQ_IMB =    models.CharField(blank=True, null=True, max_length=15, verbose_name = "IMB")	
    Vendor =    models.CharField(blank=True, null=True, max_length=15, verbose_name = "Vendor")	
    Vendor_CatNo = models.CharField(blank=True, null=True, max_length=15, verbose_name = "CatNo")	
    Salt_Form = models.CharField(blank=True, null=True, max_length=15, verbose_name = "SaltForm")	
    MW = models.DecimalField(blank=True, null=True, max_digits=12, decimal_places=2, verbose_name = "MW")	
    MW_Full = models.DecimalField(blank=True, null=True, max_digits=12, decimal_places=2, verbose_name = "MW Full")
    MF = models.CharField(blank=True, null=True, max_length=15, verbose_name = "MF")	
    MF_Full = models.CharField(blank=True, null=True, max_length=15, verbose_name = "MF Full")	
    SMILES = models.CharField(blank=True, null=True, max_length=2048, verbose_name = "SMILES")	
    smol = models.MolField(blank=True, null=True, verbose_name = "MOL")	
    ffp2 = models.BfpField(null=True)

    def __str__(self) -> str:
        return f"{self.Drug_Name} ({self.Drug_Code})"
    
    class Meta:
        db_table='drug'

#-------------------------------------------------------------------------------------------------
class VITEK_Card(AuditModel):
    """
    List of VITEK Cards
    """
#-------------------------------------------------------------------------------------------------
    Choice_Dictionary= {
        'Card_Type':'Card_Type',
    }

    Card_Barcode = models.CharField(primary_key=True, max_length=25, verbose_name = "Card Barcode") 
    OrgBatch_ID = models.ForeignKey(Organism_Batch, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING) 
    Card_Type = models.CharField(blank=False, null=False, max_length=5, verbose_name = "Card Type")
    Card_Code = models.CharField(blank=False, null=False, max_length=15, verbose_name = "Card Code") 
    Expiry_Date = models.DateField(null=True, blank=True, verbose_name = "Expiry Date")
    Instrument = models.CharField(blank=True, null=True, unique=True, max_length=50, verbose_name = "Vitek SN") 
    Proc_Date = models.DateField(null=False, blank=False, verbose_name = "Processing Date")
    Analysis_Time = models.CharField(blank=True, null=True, unique=True, max_length=15, verbose_name = "Analysis Time") 

    class Meta:
        ordering=['Card_Type','Card_Code','Card_Barcode']
        indexes = [
            models.Index(name="vcard_ctype_idx",fields=['Card_Type']),
            models.Index(name="vcard_ccode_idx",fields=['Card_Code']),
            models.Index(name="vcard_bcode_idx",fields=['Card_Barcode']),
            models.Index(name="vcard_orgbatch_idx",fields=['OrgBatch_ID']),
        ]
        db_table='VITEK_card'
    def __str__(self) -> str:
        return f"{self.OrgBatch_ID} {self.Stock_Type} {self.N_Left}"

#-------------------------------------------------------------------------------------------------
class VITEK_AST(AuditModel):
    """
     Antimicrobial Susceptibility Testing (AST) data from VITEK Cards
    """
#-------------------------------------------------------------------------------------------------

    Card_Barcode = models.ForeignKey(VITEK_Card, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING) 
    #OrgBatch_ID = models.ForeignKey(Organism_Batch, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING) 
    #Card_Code = models.CharField(blank=False, null=False, max_length=15, verbose_name = "Card Code") 
    Process = models.CharField(blank=True, null=True, unique=True, max_length=50, verbose_name = "Vitek Process")

    Drug_ID = models.ForeignKey(Drug, verbose_name = "Drug ID", on_delete=models.DO_NOTHING)
    MIC = models.CharField(null=True, blank=True, max_length=50, verbose_name = "MIC")
    BP_Profile = models.CharField(null=True, blank=True, max_length=5, verbose_name = "Break Point")
    BP_Comment = models.CharField(null=True, blank=True, max_length=120, verbose_name = "Comment")
    Selection = models.CharField(null=True, blank=True, max_length=20, verbose_name = "Selection")
    Organism = models.CharField(null=True, blank=True, max_length=120, verbose_name = "Organism")
    Filename = models.CharField(blank=True, null=True, max_length=120, verbose_name = "PDF Filename")
    Page_No = models.IntegerField(blank=True, null=True, verbose_name = "PDF PageNo")  

    class Meta:
        # ordering=['OrgBatch']
        indexes = [
            models.Index(name="vast_drugid_idx",fields=['Drug_ID']),
            models.Index(name="vast_mic_idx",fields=['MIC']),
            models.Index(name="vast_bprofile_idx",fields=['BP_Profile']),
            models.Index(name="vast_org_idx",fields=['Organism']),
            models.Index(name="vast_fname_idx",fields=['Filename']),
        ]
        db_table='VITEK_AST'
    def __str__(self) -> str:
        return f"{self.Card_Barcode} {self.ID_Organism}"

#-------------------------------------------------------------------------------------------------
class VITEK_ID(AuditModel):
    """
     Identification Testing (ID) data from VITEK Cards
    """
#-------------------------------------------------------------------------------------------------

    Card_Barcode = models.ForeignKey(VITEK_Card, verbose_name = "Card Barcode", on_delete=models.DO_NOTHING) 
    #OrgBatch_ID = models.ForeignKey(Organism_Batch, verbose_name = "OrgBatch ID", on_delete=models.DO_NOTHING) 
    #Card_Code = models.CharField(blank=False, null=False, max_length=15, verbose_name = "Card Code") 
    Process = models.CharField(blank=True, null=True, unique=True, max_length=50, verbose_name = "Vitek Process") 
    ID_Organism = models.CharField(null=True, blank=True, max_length=120, verbose_name = "ID Organism")
    ID_Probability = models.CharField(null=True, blank=True,  max_length=120, verbose_name = "ID Organism")
    ID_Confidence = models.CharField(null=True, blank=True,  max_length=120, verbose_name = "ID Organism")
    Filename = models.CharField(blank=True, null=True, max_length=120, verbose_name = "PDF Filename")
    Page_No = models.IntegerField(blank=True, null=True, verbose_name = "PDF PageNo")  

    class Meta:
        # ordering=['OrgBatch']
        indexes = [
            models.Index(name="vid_idorg_idx",fields=['ID_Organism']),
            models.Index(name="vid_idprob_idx",fields=['ID_Probability']),
            models.Index(name="vid_idconf_idx",fields=['ID_Confidence']),
            models.Index(name="vid_fname_idx",fields=['Filename']),
        ]
        db_table='VITEK_ID'
    def __str__(self) -> str:
        return f"{self.Card_Barcode}"
#-------------------------------------------------------------------------------------------------
