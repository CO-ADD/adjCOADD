import re
from model_utils import Choices
from sequences import Sequence

from rdkit import Chem
from django_rdkit import models

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError
import pgtrigger

#from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Chemical Structures and Samples 
#-------------------------------------------------------------------------------------------------

#=================================================================================================
class Chem_Structure(AuditModel):
    """
    List of ChemStructure 
    """
#=================================================================================================
    Choice_Dictionary = {
        'structure_type':'Structure_Type',
    }

    ID_SEQUENCE = 'ChemStructure'
    ID_PREFIX = 'CS'
    ID_PAD = 9

    structure_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Structure ID")
    structure_code = models.CharField(max_length=15, blank=True, verbose_name = "Structure Code")
    structure_name = models.CharField(max_length=50, blank=True, verbose_name = "Structure Name")
    structure_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="structure_type", related_name="%(class)s_structuretype")
    parent_structure_ids = ArrayField(models.CharField(max_length=15, null=True, blank=True), size=4, verbose_name = "Panel", 
                                      null=True, blank=True)
    atom_classes = ArrayField(models.CharField(max_length=15, null=True, blank=True), size=4, verbose_name = "Atom Classes", 
                                      null=True, blank=True)

    smol = models.MolField(verbose_name = "MOL")	
    tfp2 = models.BfpField(verbose_name = "Topological-Torsion FP")	
    ffp2 = models.BfpField(verbose_name = "Feature Morgan FP (FCFP)")
    mfp2 = models.BfpField(verbose_name = "Morgan FP (ECFP)")

    nfrag = models.PositiveSmallIntegerField(default=1, blank=True, verbose_name ="nFrag")
    charge = models.DecimalField(max_digits=10, decimal_places=2,  blank=True, verbose_name ="Charge")
    
    # Calculated by Trigger Function
    inchikey = models.CharField(max_length=50, blank=True,verbose_name ="InChiKey")
    mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    natoms = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nAtoms")
    hba = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="HBond Acc")
    hbd = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="HBond Don")
    logp = models.DecimalField(max_digits=9, decimal_places=2, default=0, blank=True, verbose_name ="logP")
    tpsa = models.DecimalField(max_digits=12, decimal_places=2, default=0, blank=True, verbose_name ="tPSA")
    fractioncsp3 = models.DecimalField(max_digits=9, decimal_places=2, default=0, blank=True, verbose_name ="Sp3")
    nrotbonds = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nRotBond")
    nrings = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nRings")
    narorings = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nAroRings")
    nhetarorings = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nHetAroRings")
    nhetaliphrings = models.PositiveSmallIntegerField(default=0, blank=True, verbose_name ="nHetAliphRings")

    class Meta:
        app_label = 'dchem'
        db_table = 'chem_structure'
        ordering=['structure_name']
        indexes = [
            models.Index(name="cstruct_dname_idx", fields=['structure_name']),
            models.Index(name="cstruct_dcode_idx", fields=['structure_code']),
            models.Index(name="cstruct_inchi_idx", fields=['inchikey']),
            models.Index(name="cstruct_mf_idx", fields=['mf']),
            models.Index(name="cstruct_mw_idx", fields=['mw']),
            models.Index(name="cstruct_natoms_idx", fields=['natoms']),
            models.Index(name="cstruct_nfrag_idx", fields=['nfrag']),
            models.Index(name="cstruct_charge_idx", fields=['charge']),
            GistIndex(name="cstruct_smol_idx",fields=['smol']),
            GistIndex(name="cstruct_ffp2_idx",fields=['ffp2']),
            GistIndex(name="cstruct_mfp2_idx",fields=['mfp2']),
            GistIndex(name="cstruct_tfp2_idx",fields=['tfp2'])
        ]
        triggers = [pgtrigger.Trigger(
                        name= "trigfunc_chemstruct_biu",
                        operation = pgtrigger.Insert | pgtrigger.Update,
                        when = pgtrigger.Before,
                        func = """
                                New.mfp2 := morganbv_fp(NEW.sMol);
                                New.ffp2 := featmorganbv_fp(NEW.sMol);
                                New.tfp2 := torsionbv_fp(NEW.sMol);
                                New.inchikey := mol_inchikey(NEW.sMol);
                                New.mw := mol_amw(NEW.sMol);
                                New.mf := mol_formula(NEW.sMol);
                                New.natoms := mol_numheavyatoms(NEW.sMol);
                                New.logp := mol_logp(NEW.sMol);
                                New.tpsa := mol_tpsa(NEW.sMol);
                                New.nrotbonds = mol_numrotatablebonds(NEW.sMol);
                                New.fractioncsp3 = mol_fractioncsp3(NEW.sMol);
                                New.hba = mol_hba(NEW.sMol);
                                New.hbd = mol_hbd(NEW.sMol);
                                New.nrings = mol_numrings(NEW.sMol);
                                New.narorings = mol_numaromaticrings(NEW.sMol);
                                New.nhetarorings = mol_numaromaticheterocycles(NEW.sMol);
                                New.nhetaliphrings = mol_numaliphaticheterocycles(NEW.sMol);
                                RETURN NEW;
                            """
                            )
                    ]

    #------------------------------------------------
    # def __str__(self) -> str:
    #     return f"{self.drug_id}"

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.structure_name} ({self.structure_id})"

    #------------------------------------------------
    @classmethod
    def get(cls,StructureID=None,StructureName=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if StructureID:
                retInstance = cls.objects.get(structure_id=StructureID)
            elif StructureName:
                retInstance = cls.objects.get(structure_name=StructureName)
            else:
                retInstance = None
        except:
            retInstance = None
            if verbose:
                if StructureID:
                    print(f"[Structure Not Found] {StructureID} ")
                elif StructureName:
                    print(f"[Structure Not Found] {StructureName} ")
        return(retInstance)

    @classmethod
    def get_bySmiles(cls,Smiles,verbose=0):
    # Returns an instance by smiles exact search
        try:
            retInstance = cls.objects.filter(smol__exact=Smiles).first()
        except:
            retInstance = None
            if verbose:
                print(f"[Structure Not Found] {Smiles} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,StructureID=None,StructureName=None,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        if StructureID:
            retValue = cls.objects.filter(structure_id=StructureID).exists()
        elif StructureName:
            retValue = cls.objects.filter(structure_name=StructureName).exists()
        else:
            retValue = False
        return(retValue)

    #------------------------------------------------
    @classmethod
    def exists_bySmiles(cls,Smiles,verbose=0):
    # Returns if an instance exists by drug_name or durg_id
        retValue = cls.objects.filter(smol__exact=Smiles).exists()
        return(retValue)

    #------------------------------------------------
    @classmethod
    def smiles2mol(cls,Smiles,verbose=0):
        try:
            xmol = Chem.MolFromSmiles(Smiles)
        except:
            xmol = None
            if verbose:
                print(f"[Invalid SMILES] {Smiles} ")
        return(xmol)
    
    #------------------------------------------------
    def set_molecule(self,Smiles):
        self.smol = self.smiles2mol(Smiles,verbose=0) 

    #------------------------------------------------
    def set_properties(self):
        self.charge = Chem.GetFormalCharge(self.smol)
    
    #------------------------------------------------
    def get_smiles(self):
        return(Chem.MolToSmiles(self.smol))

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if self.smol:
            self.set_properties()
            if not self.structure_id:
                self.structure_id = self.next_id()
                if self.structure_id: 
                    super(Chem_Structure, self).save(*args, **kwargs)
            else:
                super(Chem_Structure, self).save(*args, **kwargs)
        else: 
            print(f"[Not a valid Molecule] ")

# #=================================================================================================
# class Library(AuditModel):
#     """
#     List of Sample Library
#     """
# #=================================================================================================
#     Choice_Dictionary = {
#         'library_class':'Library_Class',
#     }
#     ID_SEQUENCE = 'ChemLibrary'
#     ID_PREFIX = 'LIB'
#     ID_PAD = 5

#     library_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Library ID")
#     library_name = models.CharField(max_length=50, unique=True, verbose_name = "Library Name")
#     chem_class = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Class", on_delete=models.DO_NOTHING,
#         db_column="chem_class", related_name="%(class)s_chemclass")
#     source = models.CharField(max_length=250, blank=True, verbose_name = "Source")
#     source_code = models.CharField(max_length=120, blank=True, verbose_name = "Source Code")
#     reference = models.CharField(max_length=150, blank=True, verbose_name = "Reference")

#     class Meta:
#         app_label = 'dchem'
#         db_table = 'library'
#         ordering=['library_name']
#         indexes = [
#             models.Index(name="lib_lname_idx", fields=['library_name']),
#         ]

#     #------------------------------------------------
#     def __repr__(self) -> str:
#         return f"{self.library_name} ({self.library_id}) {self.source}"

#     #------------------------------------------------
#     @classmethod
#     def get(cls,LibraryID=None,LibraryName=None,verbose=0):
#     # Returns an instance by structure_id or structure_name
#         try:
#             if LibraryID:
#                 retInstance = cls.objects.get(library_id=LibraryID)
#             elif LibraryName:
#                 retInstance = cls.objects.get(library_name=LibraryName)
#             else:
#                 retInstance = None
#         except:
#             retInstance = None
#             if verbose:
#                 if LibraryID:
#                     print(f"[Library Not Found] {LibraryID} ")
#                 elif LibraryName:
#                     print(f"[Library Not Found] {LibraryName} ")
#         return(retInstance)

#     #------------------------------------------------
#     @classmethod
#     def exists(cls,LibraryID=None,LibraryName=None,verbose=0):
#     # Returns if an instance exists by drug_name or durg_id
#         if LibraryID:
#             retValue = cls.objects.filter(library_id=LibraryID).exists()
#         elif LibraryName:
#             retValue = cls.objects.filter(library_name=LibraryName).exists()
#         else:
#             retValue = False
#         return(retValue)


#     #------------------------------------------------
#     def save(self, *args, **kwargs):
#         if not self.library_id:
#             self.library_id = self.next_id()
#             if self.library_id: 
#                 super(Library, self).save(*args, **kwargs)
#         else:
#             super(Library, self).save(*args, **kwargs) 



#=================================================================================================
class Chem_Salt(AuditModel):
    """
    List of Salt/Ion/Solvent 
    """
#=================================================================================================
    Choice_Dictionary = {
        'salt_type':'Salt_Type',
    }

    salt_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Salt ID")
    salt_code = models.CharField(max_length=15, blank=True, verbose_name = "Salt Code")
    salt_name = models.CharField(max_length=50, blank=True, verbose_name = "Salt Name")
    salt_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="salt_type", related_name="%(class)s_salttype")
    smiles = models.CharField(max_length=500, blank=True, verbose_name = "Smiles")
    smol = models.MolField(verbose_name = "MOL")	
    mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    natoms = models.IntegerField(default=0, blank=True, verbose_name ="nAtoms")
    charge = models.DecimalField(max_digits=10, decimal_places=2,  blank=True, verbose_name ="Charge")
    h_equiv = models.IntegerField(default=0, blank=True, verbose_name ="H Equivalent")

    class Meta:
        app_label = 'dchem'
        db_table = 'chem_salt'
        ordering=['salt_code']
        indexes = [
            models.Index(name="csalt_dcode_idx", fields=['salt_code']),
            models.Index(name="csalt_type_idx", fields=['salt_type']),
            GistIndex(name="csalt_smol_idx",fields=['smol']),
            # GistIndex(name="cstruct_ffp2_idx",fields=['ffp2']),
            # GistIndex(name="cstruct_mfp2_idx",fields=['mfp2']),
            # GistIndex(name="cstruct_tfp2_idx",fields=['tfp2'])
        ]
        triggers = [pgtrigger.Trigger(
                        name= "trigfunc_chemsalt_biu",
                        operation = pgtrigger.Insert | pgtrigger.Update,
                        when = pgtrigger.Before,
                        func = """
                                New.mw := mol_amw(NEW.sMol);
                                New.mf := mol_formula(NEW.sMol);
                                New.natoms := mol_numheavyatoms(NEW.sMol);
                                RETURN NEW;
                            """
                            )
                    ]

#=================================================================================================
class Chem_Alert(AuditModel):
    """
    List of Chemical Reaction/Transformations/Substructures/Alerts
    """
#=================================================================================================
    Choice_Dictionary = {
        'alert_type':'Alert_Type',
    }

    alert_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Alert ID")
    alert_code = models.CharField(max_length=15, blank=True, verbose_name = "Alert Code")
    alert_name = models.CharField(max_length=50, blank=True, verbose_name = "Alert Name")
    alert_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="alert_type", related_name="%(class)s_alerttype")
    smarts = models.CharField(max_length=10125, blank=True, verbose_name = "Smarts")
    allow_min = models.IntegerField(default=0, blank=True, verbose_name ="nAtoms")
    allow_ax = models.IntegerField(default=0, blank=True, verbose_name ="nAtoms")


    class Meta:
        app_label = 'dchem'
        db_table = 'chem_alert'
        ordering=['alert_code']
        indexes = [
            models.Index(name="calert_dcode_idx", fields=['alert_code']),
            models.Index(name="calert_type_idx", fields=['alert_type']),
            #GistIndex(name="calert_smol_idx",fields=['smol']),
            # GistIndex(name="cstruct_ffp2_idx",fields=['ffp2']),
            # GistIndex(name="cstruct_mfp2_idx",fields=['mfp2']),
            # GistIndex(name="cstruct_tfp2_idx",fields=['tfp2'])
        ]



#=================================================================================================
class Chem_Reaction(AuditModel):
    """
    List of Chemical Reaction/Transformations/Substructures/Alerts
    """
#=================================================================================================
    Choice_Dictionary = {
        'alert_type':'Alert_Type',
    }

    reaction_id = models.CharField(max_length=15, primary_key=True, verbose_name = "Reaction ID")
    reaction_code = models.CharField(max_length=15, blank=True, verbose_name = "Reaction Code")
    reaction_name = models.CharField(max_length=50, blank=True, verbose_name = "Reaction Name")
    reaction_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="reaction_type", related_name="%(class)s_reactiontype")
    smarts = models.CharField(max_length=10125, blank=True, verbose_name = "Smarts")
    run_status = models.IntegerField(default=0, blank=True, verbose_name ="RunStatus")


    class Meta:
        app_label = 'dchem'
        db_table = 'chem_reaction'
        ordering=['reaction_code']
        indexes = [
            models.Index(name="creact_dcode_idx", fields=['reaction_code']),
            models.Index(name="creact_type_idx", fields=['reaction_type']),
            #GistIndex(name="calert_smol_idx",fields=['smol']),
            # GistIndex(name="cstruct_ffp2_idx",fields=['ffp2']),
            # GistIndex(name="cstruct_mfp2_idx",fields=['mfp2']),
            # GistIndex(name="cstruct_tfp2_idx",fields=['tfp2'])
        ]


# #=================================================================================================
# class Sample(AuditModel):
#     """
#     List of ChemStructure 
#     """
# #=================================================================================================
#     ID_SEQUENCE = 'Sample'
#     ID_PREFIX = 'S'
#     ID_PAD = 9

#     sample_id = models.CharField(max_length=15,primary_key=True, verbose_name = "Sample ID")
#     sample_name = models.CharField(max_length=50, unique=True, verbose_name = "Sample Name")
#     library_id = models.ForeignKey(Library, null=False, blank=False, verbose_name = "Library ID", on_delete=models.DO_NOTHING,
#         db_column="library_id", related_name="%(class)s_library_id")

#     sample_code = models.CharField(max_length=50, unique=True, verbose_name = "Sample Code")
#     structure_id = models.ForeignKey(Chem_Structure, null=False, blank=False, verbose_name = "Structure ID", on_delete=models.DO_NOTHING,
#         db_column="structure_id", related_name="%(class)s_structure_id")
#     other_ids = models.CharField(max_length=50, unique=True, verbose_name = "Other IDs")

#     chemist = models.ForeignKey(ApplicationUser, null=True, blank=True, verbose_name = "Chemist", on_delete=models.DO_NOTHING, 
#         db_column="chemist", related_name="%(class)s_chemist")

#     class Meta:
#         app_label = 'dchem'
#         db_table = 'sample'
#         ordering=['library_id','sample_name']
#         indexes = [
#             models.Index(name="smp_sname_idx", fields=['sample_name']),
#             models.Index(name="smp_scode_idx", fields=['sample_code']),
#         ]


