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
        "drug_id":"Drug ID",
        "drug_name":"Drug Name",
        "drug_othernames":"Other Name",
        "drug_codes":"Drug Code",
        "drug_type":"Drug Type",
        "drug_note":"Drug Note"
    }

    Choice_Dictionary = {
        'drug_type':'Drug_Type',
        'max_phase':'Max_Phase',
    }

    drug_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Drug ID")
    drug_name = models.CharField(max_length=50, unique=True, verbose_name = "Drug Name")
    drug_othernames = ArrayField(models.CharField(max_length=60, blank=True),size=30, null=True)
    drug_codes = ArrayField(models.CharField(max_length=10, blank=True),size=30, null=True)
    drug_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Drug Type", on_delete=models.DO_NOTHING,
        db_column="drug_type", related_name="%(class)s_drugtype")
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
    torsionbv = models.BfpField(null=True)	
    ffp2 = models.BfpField(null=True, verbose_name = "FFP2")
    mfp2 = models.BfpField(null=True, verbose_name = "MFP2")
    #salt_form = models.CharField(blank=True, max_length=15, verbose_name = "SaltForm")	

    #------------------------------------------------
    class Meta:
        app_label = 'ddrug'
        db_table = 'drug'
        ordering=['drug_name']
        indexes = [
            models.Index(name="drug_dname_idx", fields=['drug_name']),
            GistIndex(name="drug_smol_idx",fields=['smol']),
            GistIndex(name="drug_ffp2_idx",fields=['ffp2']),
            GistIndex(name="drug_mfp2_idx",fields=['mfp2'])
        ]

    #------------------------------------------------
    def __str__(self) -> str:
        return f"{self.drug_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.drug_name} ({self.drug_id})"

   #------------------------------------------------
    @classmethod
    def str_DrugID(cls,DrugNo) -> str:
        return(f"AMD{DrugNo:05d}")

    #------------------------------------------------
    @classmethod
    def find_Next_DrugID(cls) -> str:
        Drug_IDSq=Sequence('Drug')
        Drug_nextID = next(Drug_IDSq)
        Drug_strID = cls.str_DrugID(Drug_nextID)
        while cls.exists(None,Drug_strID):
            Drug_nextID = next(Drug_IDSq)
            Drug_strID = cls.str_DrugID(Drug_nextID)
        return(Drug_strID)    

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
        cls.objects.update(ffp2=FEATMORGANBV_FP('smol'),
                           mfp2=MORGANBV_FP('smol'), 
                           torsionbv=TORSIONBV_FP('smol')
                           )

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
            self.drug_id = self.find_Next_DrugID()
            if self.drug_id: 
                super(Drug, self).save(*args, **kwargs)
                # self.__dict__.update(ffp2=FEATMORGANBV_FP('smol'), mfp2=MORGANBV_FP('smol'), torsionbv=TORSIONBV_FP('smol'))
                Drug.objects.filter(drug_id=self.drug_id).update(
                    ffp2=FEATMORGANBV_FP('smol'), 
                    mfp2=MORGANBV_FP('smol'), 
                    torsionbv=TORSIONBV_FP('smol')
                    )
        else:
            self.__dict__.update(
                ffp2=FEATMORGANBV_FP('smol'), 
                mfp2=MORGANBV_FP('smol'), 
                torsionbv=TORSIONBV_FP('smol')
                )
            super(Drug, self).save(*args, **kwargs) 
        
#=================================================================================================
class VITEK_Card(AuditModel):
#     """
#     List of VITEK Cards
#     """
#=================================================================================================
    HEADER_FIELDS = {
        "orgbatch_id":"Orgbatch",
        #"card_barcode":"Barcode",
        "card_type":"Card Type",
        "card_code":"Card Code",
        "expiry_date":"Expiry",
        "instrument":"Instrument",
        "proc_date":"Procesed",
        "analysis_time":"Time",
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
        return f"{self.card_code} ({self.card_barcode}) {self.orgbatch_id}  "
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

#    #------------------------------------------------
#     @classmethod
#     def check_from_dict(cls,cDict,valLog):
#     #
#     # Returns an instance from dictionary 
#     #  with Validation_Log for validation check
#     #  .validStatus if validated 
#     #
#         validStatus = True
       
#         retInstance = cls.get(cDict['CARD_BARCODE'])
#         if retInstance is None:
#             retInstance = cls()
#             retInstance.card_barcode = cDict['CARD_BARCODE']
#             valLog.add_log('Info','New VITEK card',f"{cDict['CARD_BARCODE']}-{cDict['CARD_CODE']}",'-')
#         else:
#             valLog.add_log('Info','Update VITEK card',f"{retInstance} -{cDict['CARD_CODE']}",'-')

#         OrgBatch = Organism_Batch.get(cDict['ORGBATCH_ID']) 
#         if OrgBatch is None:
#             valLog.add_log('Error','Organism Batch does not Exists',cDict['ORGBATCH_ID'],'Use existing OrganismBatch ID')
#             validStatus = False
#         retInstance.orgbatch_id = OrgBatch

#         retInstance.card_type = Dictionary.get(retInstance.Choice_Dictionary["card_type"],cDict['CARD_TYPE'])
#         if retInstance.card_type is None:
#             valLog.add_log('Error','Vitek Card Type not Correct',cDict['CARD_TYPE'],'-')
#             validStatus = False

#         retInstance.card_code = cDict['CARD_CODE']
#         retInstance.instrument = cDict['INSTRUMENT']
#         retInstance.expiry_date = cDict['EXPIRY_DATE']
#         retInstance.proc_date = cDict['PROCESSING_DATE']
#         retInstance.analysis_time = cDict['ANALYSIS_TIME']

#         retInstance.clean_Fields()
#         validDict = retInstance.validate()
#         if validDict:
#             validStatus = False
#             for k in validDict:
#                 valLog.add_log('Warning',validDict[k],k,'-')
#         retInstance.VALID_STATUS = validStatus
#         return(retInstance)

#=================================================================================================
class VITEK_AST(AuditModel):
    """
      Antimicrobial Suceptibility Testing (AST) data from VITEK Cards
    """
#=================================================================================================
    HEADER_FIELDS = {
        "card_barcode":"Barcode",
        "process":"Process",
        "id_organism":"ID organism",
        "id_probability":"ID Probability",
        "id_confidence":"ID Confidence",
        "id_source":"Source",
        "filename":"PDF Name",
        "page_no":"PDF PageNo"
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
    organism = models.CharField(max_length=120, blank=True, verbose_name = "Organism")
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

    #------------------------------------------------
    # @classmethod
    # def check_from_dict(cls,cDict,valLog):
    #
    # Returns an instance from dictionary 
    #  with Validation_Log for validation check
    #  .validStatus if validated 
    #
        # validStatus = True
        # Barcode = VITEK_Card.get(cDict['CARD_BARCODE']) 
        # if Barcode is None:
        #     validStatus = False
        #     valLog.add_log('Error','VITEK card does not Exists',f"{cDict['CARD_CODE']} ({cDict['CARD_BARCODE']})",'-')

        # DrugID = Drug.get(cDict['DRUG_NAME'])
        # if DrugID is None:
        #     validStatus = False
        #     valLog.add_log('Error','Drug does not Exists',f"{cDict['DRUG_NAME']} ({cDict['CARD_BARCODE']})",'-')

        # if validStatus:
        #     retInstance = cls.get(Barcode,DrugID,cDict['BP_SOURCE'],cDict['SELECTED_ORGANISM'])
        # else:
        #     retInstance = None
               
        # if retInstance is None:
        #     retInstance = cls()
        #     retInstance.card_barcode = Barcode
        #     retInstance.drug_id = DrugID
        #     retInstance.bp_source = cDict['BP_SOURCE']
        #     valLog.add_log('Info','New VITEK AST',f"{Barcode} {DrugID} {cDict['BP_SOURCE']}",'-')
        
        # retInstance.mic = cDict['MIC']
        # retInstance.process = cDict['VITEK_PROCESS']
        # retInstance.bp_profile = cDict['BP_PROFILE']
        # retInstance.bp_comment = cDict['BP_COMMENT']
        # retInstance.selection = cDict['ORGANISM_ORIGIN']
        # retInstance.organism = cDict['SELECTED_ORGANISM']
        # retInstance.filename = cDict['FILENAME']
        # retInstance.page_no = cDict['PAGENO']  

        # retInstance.clean_Fields()
        # validDict = retInstance.validate()
        # if validDict:
        #     validStatus = False
        #     for k in validDict:
        #         valLog.add_log('Warning',validDict[k],k,'-')

        # retInstance.VALID_STATUS = validStatus
        # return(retInstance)
    
#=================================================================================================
class VITEK_ID(AuditModel):
    """
      Identification Testing (ID) data from VITEK Cards
    """
#=================================================================================================
    HEADER_FIELDS = {
        "card_barcode":"Barcode",
        "drug_id":"Drug",
        "mic":"MIC",
        "process":"Vitek Process",
        "bp_profile":"Break Point",
        "bp_comment":"Comment",
        "bp_source":"Source",
        "selection":"Selection",
        "organism":"Organism",
        "filename":"PDF Filename",
        "page_no":"PDF pageNo"
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
        #ordering=['card_barcode']
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

   #------------------------------------------------
    # @classmethod
    # def check_from_dict(cls,cDict,valLog):
    # #
    # # Returns an instance from dictionary 
    # #  with Validation_Log for validation check
    # #  .validStatus if validated 
    # #
    #     validStatus = True
    #     Barcode = VITEK_Card.get(cDict['CARD_BARCODE']) 
    #     if Barcode is None:
    #         validStatus = False
    #         valLog.add_log('Error','VITEK card does not Exists',f"{cDict['CARD_CODE']} ({cDict['CARD_BARCODE']})",'-')

    #     retInstance = cls.get(Barcode)
    #     if retInstance is None:
    #         retInstance = cls()
    #         retInstance.card_barcode = Barcode
    #         valLog.add_log('Info','New VITEK ID',f"{cDict['CARD_CODE']} ({cDict['CARD_BARCODE']})",'-')
    #     else:
    #         valLog.add_log('Info','Update VITEK ID',f"{cDict['CARD_CODE']} ({Barcode})",'-')

    #     retInstance.process = cDict['VITEK_PROCESS']
    #     retInstance.id_organism = cDict['ID_ORGANISM']
    #     retInstance.id_probability = cDict['ID_PROBABILITY']
    #     retInstance.id_confidence = cDict['ID_CONFIDENCE']
    #     #retInstance.id_source = cDict['CARD_BARCODE']
    #     retInstance.filename = cDict['FILENAME']
    #     retInstance.page_no = cDict['PAGENO']  

    #     retInstance.clean_Fields()
    #     validDict = retInstance.validate()
    #     if validDict:
    #         validStatus = False
    #         for k in validDict:
    #             valLog.add_log('Warning',validDict[k],k,'-')

    #     retInstance.VALID_STATUS = validStatus
    #     return(retInstance)

#=================================================================================================
class MIC_COADD(AuditModel):
    """
     Antibiogram from CO-ADD screening    
    """
#=================================================================================================
    HEADER_FIELDS = {
        # example fields for test view
        "drug_id.drug_name":"Drug Name",
        "mic_type":"Type",
        "mic":"MIC",
        "orgbatch_id.organism_id.gen_property":"Organism Resistance Property",
        "run_id":"Run ID",
        "bp_profile":"Break Point",
        "media":"Media",
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

    # Future update to ForeignKey (JZG) ... same for Testplate_ID
    #run_id = models.ForeignKey(Screen_Run, null=True, blank=True, verbose_name = "Run ID", on_delete=models.DO_NOTHING,
    #    db_column="run_id", related_name="%(class)s_runid")
    run_id = models.CharField(max_length=25, blank=True, verbose_name = "RunID")
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
        retStr = f"{self.drug_id.drug_name} {self.orgbatch_id} {self.mic} {self.run_id}"
        return(retStr)

   #------------------------------------------------
    @classmethod
    def get(cls,OrgBatchID,DrugID,TestPlateID,TestWellID,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        try:
            retInstance = cls.objects.get(orgbatch_id=OrgBatchID,run_id=DrugID,testplate_id=TestPlateID,testwell_id=TestWellID)
        except:
            if verbose:
                print(f"[MIC Not Found] {OrgBatchID} {DrugID} {TestPlateID} {TestWellID}")
            retInstance = None
        return(retInstance)

   #------------------------------------------------
    @classmethod
    def exists(cls,OrgBatchID,DrugID,TestPlateID,TestWellID,verbose=0):
    # Returns an instance if found by OrgBatchID and DrugID
        return cls.objects.filter(orgbatch_id=OrgBatchID,run_id=DrugID,testplate_id=TestPlateID,testwell_id=TestWellID).exists()

    # #------------------------------------------------
    # @classmethod
    # def check_from_dict(cls,cDict,valLog):
    # #
    # # Returns an instance from dictionary 
    # #  with Validation_Log for validation check
    # #  .validStatus if validated 
    # #
    #     validStatus = True
    #     DrugID = Drug.get(cDict['DRUG_NAME'])
    #     if DrugID is None:
    #         validStatus = False
    #         valLog.add_log('Error','Drug does not Exists',f"{cDict['DRUG_NAME']} ",'-')

    #     OrgBatchID = Organism_Batch.get(cDict['ORGBATCH_ID']) 
    #     if OrgBatchID is None:
    #         validStatus = False
    #         valLog.add_log('Error','OrgBatchID does not Exists',f"{cDict['ORGBATCH_ID']} ",'-')

    #     if validStatus:
    #         retInstance = cls.get(OrgBatchID,DrugID,cDict['TESTPLATE_ID'],cDict['TESTWELL_ID'])
    #     else:
    #         retInstance = None
               
    #     if retInstance is None:
    #         retInstance = cls()
    #         retInstance.orgbatch_id = OrgBatchID
    #         retInstance.drug_id = DrugID
    #         retInstance.run_id = cDict['RUN_ID']
    #         retInstance.testplate_id = cDict['TESTPLATE_ID']
    #         retInstance.testwell_id = cDict['TESTWELL_ID']
    #         valLog.add_log('Info','New MIC ',f"{OrgBatchID} {DrugID} {cDict['TESTPLATE_ID']}:{cDict['TESTWELL_ID']}",'-')
        
    #     retInstance.mic = cDict['MIC']
    #     retInstance.mic_unit = cDict['MIC_UNIT']
    #     retInstance.mic_type = Dictionary.get(cls.Choice_Dictionary["mic_type"],'BMD',None,verbose=1)

    #     retInstance.plate_size = Dictionary.get(cls.Choice_Dictionary["plate_size"],cDict['PLATE_SIZE'],None,verbose=1)
    #     retInstance.plate_material = Dictionary.get(cls.Choice_Dictionary["plate_material"],cDict['PLATE_MATERIAL'],None,verbose=1)

    #     #retInstance.bp_profile = cDict['BP_PROFILE']
    #     #retInstance.bp_source = cDict['BP_SOURCE']
    #     #retInstance.media = Dictionary.get(cls.Choice_Dictionary["media"],cDict['MEDIA'],None,verbose=1)

    #     retInstance.clean_Fields()
    #     validDict = retInstance.validate()
    #     if validDict:
    #         validStatus = False
    #         for k in validDict:
    #             valLog.add_log('Warning',validDict[k],k,'-')

    #     retInstance.VALID_STATUS = validStatus
        
    #     return(retInstance)
    
    # overide get value and get value from parent model
    def get_values(self, fields=None):
        from django.db.models import Model
        if fields is None:
            fields = self.HEADER_FIELDS

        value_list=[]
        fieldsname=[field.name for field in self._meta.fields]
        for name in fields.keys():
            n=len(name.split("."))
            nameArray=name.split(".")
            if n>1 and nameArray[0] in fieldsname:
                obj_parent=getattr(self, nameArray[0])
                obj = obj_parent
                i=1
                while i < n:
                    obj = getattr(obj, name.split(".")[i])
                    i += 1
                    if i>=5:
                        
                        break
                value_list.append(obj)
                
            elif name in fieldsname:
                obj=getattr(self, name)
                if obj:
                    if isinstance(obj, Model):
                        value_list.append(obj.pk)
                    elif isinstance(obj, list):
                        array_to_string=','.join(str(e) for e in obj)
                        value_list.append(array_to_string) 
                    else:   
                        value_list.append(obj)
                else:
                    value_list.append(" ")
                    
        return value_list

    # override to get field contains foreignkey fields
    @classmethod
    def get_fields(cls, fields=None):
        if fields is None:
            fields = cls.HEADER_FIELDS
        if fields:
            fieldsname=[field.name for field in cls._meta.fields]
            select_fields=[fields[f] for f in fields.keys() if f in fieldsname or f.split(".")[0] in fieldsname]
        else:
            select_fields=None   
        return select_fields

#=================================================================================================
class MIC_Pub(AuditModel):
    """
     Antibiogram from Public sources    
    """
#=================================================================================================
    HEADER_FIELDS   = {
           # example fields for test view
        "drug_id.drug_name":"Drug Name",
        "mic_type":"Type",
        "mic":"MIC",
        "orgbatch_id.organism_id.gen_property":"Organism Resistance Property",
    
        "source":"Source",
        "bp_profile":"Break Point",
    }
    Choice_Dictionary = {
        'mic_type':'MIC_Type',
    }

    organism_id = models.ForeignKey(Organism, null=False, blank=False, verbose_name = "Organism ID", on_delete=models.DO_NOTHING,
        db_column="oragnism_id", related_name="%(class)s_organism_id") 
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

    # #------------------------------------------------
    # @classmethod
    # def check_from_dict(cls,cDict,valLog):
    # #
    # # Returns an instance from dictionary 
    # #  with Validation_Log for validation check
    # #  .validStatus if validated 
    # #
    #     #print(cDict)
    #     validStatus = True
    #     DrugID = Drug.get(cDict['DRUG_NAME'])
    #     if DrugID is None:
    #         validStatus = False
    #         valLog.add_log('Error','Drug does not Exists',f"{cDict['DRUG_NAME']} ",'-')

    #     OrganismID = Organism.get(cDict['ORGANISM_ID']) 
    #     if OrganismID is None:
    #         validStatus = False
    #         valLog.add_log('Error','OrganismID does not Exists',f"{cDict['ORGANISM_ID']} ",'-')

    #     if validStatus:
    #         retInstance = cls.get(OrganismID,DrugID,cDict['SOURCE'])
    #     else:
    #         retInstance = None
               
    #     if retInstance is None:
    #         retInstance = cls()
    #         retInstance.organism_id = OrganismID
    #         retInstance.drug_id = DrugID
    #         retInstance.source = cDict['SOURCE']
    #         valLog.add_log('Info','New MIC ',f"{OrganismID} {DrugID} {cDict['SOURCE']}",'-')
        
    #     retInstance.mic = cDict['MIC']
    #     retInstance.mic_unit = cDict['MIC_UNIT']
    #     retInstance.mic_type = Dictionary.get(cls.Choice_Dictionary["mic_type"],cDict['SOURCE_TYPE'],None,verbose=1)
    #     retInstance.bp_profile = cDict['BP_PROFILE']
    #     retInstance.bp_source = cDict['BP_SOURCE']

    #     retInstance.clean_Fields()
    #     validDict = retInstance.validate()
    #     if validDict:
    #         validStatus = False
    #         for k in validDict:
    #             valLog.add_log('Warning',validDict[k],k,'-')

    #     retInstance.VALID_STATUS = validStatus
        
    #     return(retInstance)


         # overide get value and get value from parent model
    def get_values(self, fields=None):
        from django.db.models import Model
        if fields is None:
            fields = self.HEADER_FIELDS

        value_list=[]
        fieldsname=[field.name for field in self._meta.fields]
        for name in fields.keys():
            n=len(name.split("."))
            nameArray=name.split(".")
            if n>1 and nameArray[0] in fieldsname:
                obj_parent=getattr(self, nameArray[0])
                obj = obj_parent
                i=1
                while i < n:
                    obj = getattr(obj, name.split(".")[i])
                    i += 1
                    if i>=5:
                        
                        break
                value_list.append(obj)
                
            elif name in fieldsname:
                obj=getattr(self, name)
                if obj:
                    if isinstance(obj, Model):
                        value_list.append(obj.pk)
                    elif isinstance(obj, list):
                        array_to_string=','.join(str(e) for e in obj)
                        value_list.append(array_to_string) 
                    else:   
                        value_list.append(obj)
                else:
                    value_list.append(" ")
                    
        return value_list


    # override to get field contains foreignkey fields
    @classmethod
    def get_fields(cls, fields=None):
        if fields is None:
            fields = cls.HEADER_FIELDS
        if fields:
            fieldsname=[field.name for field in cls._meta.fields]
            select_fields=[fields[f] for f in fields.keys() if f in fieldsname or f.split(".")[0] in fieldsname]
        else:
            select_fields=None   
        return select_fields