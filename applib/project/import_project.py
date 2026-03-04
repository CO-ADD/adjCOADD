import os
import pandas as pd
#from django.apps import apps
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from dcollab.models import Collab_Group, Collab_User, Organisation
from dsample.models import Project, Project_Membership,COADD_Compound
from apputil.models import Dictionary

from rdkit import Chem

#-----------------------------------------------------------------------------
def get_CompoundSubmisssion_xlsx(xlsFile, FillNA='-', **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    SHEET_NAMES = {
        'Samples': 'Sample Info',
        'Contacts': 'Collaborator Info',
    }
    
    _Sheets = {
        'Samples': None,
        'Contacts': None,
        }

    print(f"  [get_CompoundSubmisssion_xlsx] {xlsFile} ")
    if os.path.isfile(xlsFile):
        fXlsx = open(xlsFile, "rb")
        xls = pd.ExcelFile(fXlsx)
        xlsheets = xls.sheet_names
        #xls = pd.ExcelFile(xlsFile)
            
        for key in SHEET_NAMES:
            if SHEET_NAMES[key] in xlsheets:
                _Sheets[key] = xls.parse(SHEET_NAMES[key], header = None) 
            else:
                if valLog:
                    valLog.add_error(f'Missing Sheet for {key}',f"XLSX {os.path.basename(xlsFile)}",f"Correct XLSX Sheets [{SHEET_NAMES[key]}]")
        fXlsx.close()
        
    return(_Sheets)

#-----------------------------------------------------------------------------
def validate_smiles(smiles_string, *args, **kwargs):
#-----------------------------------------------------------------------------
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return mol
    except Exception:
        print(smiles_string)
        if valLog:
            valLog.add_error('Wrong SMILES',smiles_string,'Smiles not valid','Correct Smiles')
        return None

#-----------------------------------------------------------------------------
def parse_SampleInfo_Sheet(xSheet, *args, **kwargs):
#-----------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    SAMPLE_FIELDS = {
        0:'compound_code',
        1:'plate_id',
        2:'well_id',
        3:'barcode',
        4:'reg_amount',
        5:'reg_volume',
        6:'reg_conc',
        7:'reg_conc_unit',
        8:'reg_solvent',
        9:'reg_mw',
        10:'reg_mf',
        11:'reg_smiles',
        12:'structure_notes',
        13:'compound_type',
        14:'compound_note',
    }
    SAMPLE_FLOATS = ['reg_amount','reg_conc','reg_mw']
    
    _Samples = []
    _missing_smiles = 0
    _unique_codes = {}
    
    if xSheet is not None:    
        for r in range(3,len(xSheet)):
            _Sample = {}
            for i in SAMPLE_FIELDS.keys():
                _Sample[SAMPLE_FIELDS[i]] = xSheet[i][r]
                
            # Set Units
            _Sample['reg_amount_unit'] = 'mg'
            _Sample['reg_volume_unit'] = 'uL'
            
            _s = _Sample['reg_amount_unit'].split(' ')
            if len(_s)>1:
                _Sample['reg_amount'] = _s[0]
                _Sample['reg_amount_unit'] = _s[1]

            #Clean Floats
            for f in SAMPLE_FLOATS:
                try:
                    _Sample[f] = float(_Sample[f])
                except ValueError:
                    if valLog:
                        valLog.add_error("Not Numeric",_Sample[f],f" Field [{f}] in row {r} is not numeric","Correct Compound Sheet")   
            
            # Check Unique Codes
            if _Sample['compound_code'] in _unique_codes:
                valLog.add_error('Duplicate Code',_Sample['compound_code'],'Dublicate Codes in Submission','Correct Compound Code')
            else:
                _unique_codes[_Sample['compound_code']] = 1
            
            # Check Smiles
            if pd.isna(_Sample['reg_smiles']):
                _missing_smiles += 1
            else:
                _Sample['smol'] = validate_smiles(_Sample['reg_smiles'])
                if _Sample['smol'] is None:
                    if valLog:
                        valLog.add_error('Wrong SMILES',_Sample['reg_smiles'],'Smiles not valid','Correct Smiles')
                        
            # Check Types
            if pd.isna(_Sample['compound_type']):
                _Sample['compound_type'] = 'Small molecule'
            else:
                djDict = Dictionary.get_bysimilarity(COADD_Compound.DICTIONARY_FIELDS['compound_type'],_Sample['compound_type'])
                if djDict is None:
                    valLog.add_error('Unknown Compound Type',_Sample['compound_type'],'','Correct Compound Type')
                    _Sample['compound_type'] = None
                else:
                    _Sample['compound_type'] = str(djDict)

                            
            _Samples.append(_Sample)

        if _missing_smiles > 0:
            valLog.add_warning("Missing SMILES",f" {_missing_smiles} Samples without SMILES ") 
            
    valLog.add_info("Sample Info",f" {len(_Samples)} Samples") 
    return(_Samples)

#-----------------------------------------------------------------------------
def parse_ContactInfo_Sheet(xSheet, *args, **kwargs):
#-----------------------------------------------------------------------------
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    # - Project Title --------------------------------
    PRJTITLE_ROW = 3
    PRJTITLE_COL = 3

    _PrjTitle = xSheet[PRJTITLE_COL-1][PRJTITLE_ROW-1]
    if pd.isna(_PrjTitle):
        _PrjTitle = None

    # - Contacts --------------------------------
    CONTACT_FIELDS = {
        0:'title',
        1:'first_name',
        2:'last_name',
        3:'position',
        4:'email',
        5:'phone',
        6:'organisation',
        7:'department',
        8:'street_address',
        9:'country',     
    }
    CONTACT_ROW = 6
    CONTACT_COL = 3
    CONTACT_TYPES = {
        0:'LI',
        1:'PC',
        2:'AC',        
    }

    _Contacts = {}            
    for c in range(3):
        _Contact = {'type':CONTACT_TYPES[c]}
        
        # Col C(2), F(5), I(8)
        for i in CONTACT_FIELDS.keys():
            #print(f" {CONTACT_FIELDS[i]} {xSheet[(c*3)+CONTACT_COL-1][i+CONTACT_ROW-1]}")
            _Contact[CONTACT_FIELDS[i]] = xSheet[(c*3)+CONTACT_COL-1][i+CONTACT_ROW-1]
        
        # Check Email
        if pd.isna(_Contact['email']):
            if not (pd.isna(_Contact['first_name']) and pd.isna(_Contact['last_name'])):
                valLog.add_warning("Has no EMail",_Contact['email'],f" Email of [{_Contact['type']}] missing","Correct Contact Sheet") 
        else:
            try:
                validate_email(_Contact['email'])
            except ValidationError as e:
                if valLog:
                    valLog.add_error("Not EMail",_Contact['email'],f" Email of [{_Contact['type']}] not valid","Correct Contact Sheet")   

            _Contacts[CONTACT_TYPES[c]] = _Contact
        
    valLog.add_info("Collaborator Info",f" {len(_Contacts)} Contacts {list(_Contacts.keys())}")
    if _PrjTitle:
        valLog.add_info("Project Title",f" {_PrjTitle}")
        
    return(_Contacts,_PrjTitle)


#-----------------------------------------------------------------------------
def Upload_Project_Collab(djProject, CollabDict,  upload=False, overwrite=False, **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
    
    for key in CollabDict:
        
        # -- Organisation -----------------------------
        djOrg = Organisation.get_bysimilarity(CollabDict[key]['organisation'])
        
        if djOrg is None:
            djOrg = Organisation()
            djOrg.organisation_name = CollabDict[key]['organisation']
            valLog.add_warning("New Organisation",CollabDict[key]['organisation'])
            
            # - Validation ----------------------------------------------------------
            validStatus = True
            djOrg.set_defaults_model()
            validDict = djOrg.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('Warning',k,validDict[k],'-')
            if validStatus:
                if upload:
                    djOrg.save()
        else:
            valLog.add_info("Existing Organisation",f"{djOrg}")
        
        # -- Collaborator -----------------------------
        djUsr = Collab_User.get(None,CollabDict[key]['email'])
        if djUsr is None:
            djUsr = Collab_User()   
            djUsr.email = CollabDict[key]['email'] 
            djUsr.first_name = CollabDict[key]['first_name'] 
            djUsr.last_name = CollabDict[key]['last_name']
            djUsr.organisation_id = djOrg 
            valLog.add_warning("New Collaborator",f"{djUsr}")

            # - Validation ----------------------------------------------------------
            validStatus = True
            djUsr.set_defaults_model()
            validDict = djUsr.validate_fields()
            if validDict:
                validStatus = False
                for k in validDict:
                    print('Warning',k,validDict[k],'-')
            if validStatus:
                if upload:
                    djUsr.save()
        else:
            valLog.add_info("Existing Collaborator",f"{djUsr}")

        # -- Collaborator Group  -----------------------------    
        if 'LI' == key:
            djGrp = Collab_Group.get(None,Code=None, LI_ID=djUsr, Organisation_ID=djOrg)
            if djGrp is None:
                djGrp = Collab_Group()
                djGrp.group_code = f"{CollabDict[key]['last_name']}{CollabDict[key]['first_name'][0]}_{djOrg.organisation_code}"
                djGrp.organisation_id = djOrg
                valLog.add_warning("New Group",f"{djGrp}")

                # - Validation ----------------------------------------------------------
                validStatus = True
                djGrp.set_defaults_model()
                validDict = djGrp.validate_fields()
                if validDict:
                    validStatus = False
                    for k in validDict:
                        print('Warning',k,validDict[k],'-')
                if validStatus:
                    if upload:
                        djGrp.save()
            else:
                valLog.add_info("Existing Group",f"{djGrp}")
        else:
            # -- Project Membership [PC,AC] -----------------------------
            djMemb = Project_Membership.get(djUsr,djProject)
            if djMemb is None:
                djMemb = Project_Membership()
                djMemb.user_id = djUsr
                djMemb.project_id = djProject
                djMemb.role = key
            else:
                djMemb.role = key
                
            if upload:
                djMemb.save()    

            
# #-----------------------------------------------------------------------------
# def parse_CompoundSubmission(SubmissionFile, SubType=['Samples','Contacts'], **kwargs):
# #-----------------------------------------------------------------------------
#     valLog = kwargs.get('valLog',None)
#     verbose = kwargs.get('verbose',0)
    
#     _subSheets = get_CompoundSubmisssion_xlsx(SubmissionFile,Sheets=SubType,FillNA=None) 

#     # Parse Sample Submission
#     if 'Samples' in _subSheets:
#         Sample_Lst = parse_SampleSubmission(_subSheets['Samples'], valLog=valLog)
             
#     if 'Contacts' in _subSheets:
#         Contact_Lst = parse_ContactSubmission(_subSheets['Contacts'], valLog=valLog)
 

        

# def upload_CompoundSubmission_Process(Request, DirName, ExcelFile, ProjectID=None,upload=False,appuser=None):
#     #upload_CompoundSubmission
#     pass

# def upload_CompoundSubmission(ProjectID=None,upload=False,appuser=None,valLog=None):
#     #process_CompoundSubmission
#     pass

# def process_CompoundSubmission(ProjectID=None,upload=False,appuser=None,valLog=None):
#     #parse_CompoundSubmission
#     pass


# def Project_fromDict(iDict,valLog,upload=False):
#     pass

# def Compound_fromDict(iDict,valLog,upload=False):
#     pass

# def CollabUser_fromDict(iDict,valLog,upload=False):
#     pass

# def CollabGroup_fromDict(iDict,valLog,upload=False):
#     pass

# def CollabUser_fromDict(iDict,valLog,upload=False):
#     pass

# def Organisation_fromDict(iDict,valLog,upload=False):
#     pass

# def Test():
#     pass