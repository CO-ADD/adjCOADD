import re
from model_utils import Choices
from sequences import Sequence
from rdkit import Chem
from django_rdkit import models
from django_rdkit.models import MOL_TO_SMILES

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GistIndex
from django.db import transaction, IntegrityError
import pgtrigger

#from adjcoadd.constants import *
from apputil.models import AuditModel, Dictionary, ApplicationUser, Document

from adjcoadd.constants import *

#-------------------------------------------------------------------------------------------------
# Chemical Structures with Salt, Reactions and Alerts
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
    @classmethod
    def register_fromDict(cls,sDict,smiles_name='smiles',mol_name='mol',verbose=0):
        djchem = None
        if smiles_name in sDict:
            djchem = cls.get_bySmiles(sDict[smiles_name])
            if not djchem:
                djchem = cls()
                djchem.set_molecule(sDict[smiles_name])
        elif mol_name in sDict:
            djchem = cls.get_bySmiles(Chem.MolToSmiles(sDict[mol_name]))
            if not djchem:
                djchem = cls()
                djchem.smol = sDict[mol_name]

        if djchem:        
            for s in ['structure_code','structure_name']:
                if s in sDict:
                    setattr(djchem,s,sDict[s])

            validStatus = True
            djchem.clean_Fields()
            validDict = djchem.validate()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('[reg Chem_Structure] Warning',k,validDict[k])

            if validStatus:
                djchem.save()
                return(djchem)
            else:
                return(None)
        else:
            if verbose:
                print(f"[reg Chem_Structure] no {smiles_name} or {mol_name} in {sDict}")
            return(None)

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
    salt_name = models.CharField(max_length=100, blank=True, verbose_name = "Salt Name")
    salt_type = models.ForeignKey(Dictionary, null=True, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="salt_type", related_name="%(class)s_salttype")
    smiles = models.CharField(max_length=500, blank=True, verbose_name = "Smiles")
    mf = models.CharField(max_length=500, blank=True, verbose_name = "MF")
    mw = models.DecimalField(max_digits=12, decimal_places=3, default=0, blank=True, verbose_name ="MW")
    natoms = models.IntegerField(default=0, blank=True, verbose_name ="nAtoms")
    charge = models.DecimalField(max_digits=10, decimal_places=2,  blank=True, verbose_name ="Charge")
    h_equiv = models.IntegerField(default=0, blank=True, verbose_name ="H Equivalent")

    class Meta:
        app_label = 'dchem'
        db_table = 'chem_salt'
        ordering=['salt_type','salt_id']
        indexes = [
            models.Index(name="csalt_sid_idx", fields=['salt_id']),
            models.Index(name="csalt_type_idx", fields=['salt_type']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.salt_id} ({self.salt_type})"

    #------------------------------------------------
    @classmethod
    def get(cls,SaltID,SaltType=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if SaltType:
                retInstance = cls.objects.get(salt_id=SaltID,salt_type=SaltType)
            else:
                retInstance = cls.objects.get(salt_id=SaltID)
        except:
            retInstance = None
            if verbose:
                if SaltType:
                    print(f"[Salt Not Found] {SaltID} for {SaltType} ")
                else :
                    print(f"[Salt Not Found] {SaltID} ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,SaltID,verbose=0):
    # Returns if an instance exists by salt_id
        retValue = cls.objects.filter(salt_id=SaltID).exists()
        return(retValue)

    #------------------------------------------------
    @classmethod
    def get_Salt_DictList(cls,SaltType):
        return(cls.objects.filter(salt_type=SaltType).values())

    #------------------------------------------------
    @classmethod
    def get_Salt(cls,SaltID):
        return(cls.objects.filter(salt_id=SaltID).values())

 
#=================================================================================================
class Chem_Group(AuditModel):
    """
    List of Chemical Reaction/Transformations/Substructures/Alerts
    """
#=================================================================================================
    ID_SEQUENCE = 'ChemGroup'
    ID_PREFIX = 'CG'
    ID_PAD = 5

    Choice_Dictionary = {
        'chemgroup_type':'ChemGroup_Type',
    }

    chemgroup_id = models.CharField(max_length=15, primary_key=True, verbose_name = "ChemGroup ID")
    chemgroup_code = models.CharField(max_length=15, unique=True, verbose_name = "ChemGroup Code")
    chemgroup_name = models.CharField(max_length=150, blank=True, verbose_name = "ChemGroup Name")
    chemgroup_set = models.CharField(max_length=150, blank=False, verbose_name = "ChemGroup Set")
    chemgroup_type = models.ForeignKey(Dictionary, null=False, blank=True, verbose_name = "Type", on_delete=models.DO_NOTHING,
        db_column="chemgroup_type", related_name="%(class)s_chemgroup_type")
    smarts = models.CharField(max_length=10125, blank=False, verbose_name = "Smarts")
    allowed_min = models.IntegerField(default=0, blank=True, verbose_name ="Min")
    allowed_max = models.IntegerField(default=0, blank=True, verbose_name ="Max")


    class Meta:
        app_label = 'dchem'
        db_table = 'chem_group'
        ordering=['chemgroup_set','chemgroup_code']
        indexes = [
#            models.Index(name="cgrp_dcode_idx", fields=['chemgroup_code']),
            models.Index(name="cgrp_type_idx", fields=['chemgroup_type']),
            models.Index(name="cgrp_set_idx", fields=['chemgroup_set']),
            models.Index(name="cgrp_name_idx", fields=['chemgroup_name']),
        ]

    #------------------------------------------------
    def __repr__(self) -> str:
        return f"{self.chemgroup_id} {self.chemgroup_code} ({self.chemgroup_set})"

    #------------------------------------------------
    @classmethod
    def get(cls,ChemGroupID,ChemGroupCode=None,ChemGroupSet=None,verbose=0):
    # Returns an instance by structure_id or structure_name
        try:
            if ChemGroupID:
                retInstance = cls.objects.get(chemgroup_id=ChemGroupID)
            elif ChemGroupCode:
                if ChemGroupSet:
                    retInstance = cls.objects.get(chemgroup_code=ChemGroupCode,chemgroup_set=ChemGroupSet)
                else:
                    retInstance = cls.objects.get(chemgroup_code=ChemGroupCode)
        except:
            retInstance = None
            if verbose:
                if ChemGroupID:
                    print(f"[ChemGroup Not Found] {ChemGroupID} ")
                else :
                    print(f"[ChemGroup Not Found] {ChemGroupCode} ({ChemGroupSet}) ")
        return(retInstance)

    #------------------------------------------------
    @classmethod
    def exists(cls,ChemGroupID,verbose=0):
    # Returns if an instance exists by salt_id
        retValue = cls.objects.filter(chemgroup_id=ChemGroupID).exists()
        return(retValue)

    #------------------------------------------------
    def save(self, *args, **kwargs):
        if not self.chemgroup_id:
            self.chemgroup_id = self.next_id()
            if self.chemgroup_id: 
                super(Chem_Group, self).save(*args, **kwargs)
        else:
            super(Chem_Group, self).save(*args, **kwargs)

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
