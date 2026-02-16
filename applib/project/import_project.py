import os
import pandas as pd
#from django.apps import apps
from django.core.validators import validate_email
from django.core.exceptions import ValidationError

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
        4:'amount',
        5:'solution',
        6:'conc',
        7:'conc_unit',
        8:'solubility',
        9:'MW',
        10:'MF',
        11:'smiles',
        12:'structure_notes',
        13:'compound_type',
        14:'compound_note',
    }
    SAMPLE_FLOATS = ['amount','conc','MW']
    
    _Samples = []
    _missing_smiles = 0
    _unique_codes = {}
    
    if xSheet is not None:    
        for r in range(3,len(xSheet)):
            _Sample = {}
            for i in SAMPLE_FIELDS.keys():
                _Sample[SAMPLE_FIELDS[i]] = xSheet[i][r]
                
            #Clean Floats
            for f in SAMPLE_FLOATS:
                try:
                    _Sample[f] = float(_Sample[f])
                except ValueError:
                    if valLog:
                        valLog.add_error("Not Numeric",_Sample[f],f" Field [{f}] in row {r} is not numeric","Correct Compound Sheet")   
            
            # Check Unique Codes
            if _Sample['compound_code'] in _unique_codes:
                valLog.add_error('Duplicate Code',_Sample['compound_code'],'Code exists already','Correct Compound Code')
            else:
                _unique_codes[_Sample['compound_code']] = 1
            
            # Check Smiles
            if pd.isna(_Sample['smiles']):
                _missing_smiles += 1
            else:
                _Sample['smol'] = validate_smiles(_Sample['smiles'])
                if _Sample['smol'] is None:
                    if valLog:
                        valLog.add_error('Wrong SMILES',_Sample['smiles'],'Smiles not valid','Correct Smiles')
                            
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
        0:'PI',
        1:'PC',
        2:'M',        
    }

    _Contacts = {}            
    for c in range(3):
        _Contact = {'type':CONTACT_TYPES[c]}
        
        # Col C(2), F(5), I(8)
        for i in CONTACT_FIELDS.keys():
            print(f" {CONTACT_FIELDS[i]} {xSheet[(c*3)+CONTACT_COL-1][i+CONTACT_ROW-1]}")
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
def parse_CompoundSubmission(SubmissionFile, SubType=['Samples','Contacts'], **kwargs):
#-----------------------------------------------------------------------------
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
    
    _subSheets = get_CompoundSubmisssion_xlsx(SubmissionFile,Sheets=SubType,FillNA=None) 

    # Parse Sample Submission
    if 'Samples' in _subSheets:
        Sample_Lst = parse_SampleSubmission(_subSheets['Samples'], valLog=valLog)
             
    if 'Contacts' in _subSheets:
        Contact_Lst = parse_ContactSubmission(_subSheets['Contacts'], valLog=valLog)
 

        

def upload_CompoundSubmission_Process(Request, DirName, ExcelFile, ProjectID=None,upload=False,appuser=None):
    #upload_CompoundSubmission
    pass

def upload_CompoundSubmission(ProjectID=None,upload=False,appuser=None,valLog=None):
    #process_CompoundSubmission
    pass

def process_CompoundSubmission(ProjectID=None,upload=False,appuser=None,valLog=None):
    #parse_CompoundSubmission
    pass


def Project_fromDict(iDict,valLog,upload=False):
    pass

def Compound_fromDict(iDict,valLog,upload=False):
    pass

def CollabUser_fromDict(iDict,valLog,upload=False):
    pass

def CollabGroup_fromDict(iDict,valLog,upload=False):
    pass

def CollabUser_fromDict(iDict,valLog,upload=False):
    pass

def Organisation_fromDict(iDict,valLog,upload=False):
    pass

def Test():
    pass