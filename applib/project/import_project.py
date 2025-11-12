import os
import pandas as pd
#from django.apps import apps
from django.core.validators import validate_email
from django.core.exceptions import ValidationError

from rdkit import Chem

#-----------------------------------------------------------------------------
def get_CompoundSubmisssion_xlsx(xlsFile, Sheets=[], FillNA='-', **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    SHEET_NAMES = {
        'Samples': 'Compound Submission',
        'Contacts': 'Contact Info',
        
    }
    Submission_Sheets = {
        'Samples': None,
        'Contacts': None,
        }

    if os.path.isfile(xlsFile):
        fXlsx = open(xlsFile, "rb")
        xls = pd.ExcelFile(fXlsx)
        #xls = pd.ExcelFile(xlsFile)

        if Sheets is None:
            sheets = [Submission_Sheets.keys()]
        for key in Sheets:
            if SHEET_NAMES[key] in Submission_Sheets:
                Submission_Sheets[key] = xls.parse(key, header = None) 
            else:
                if valLog:
                    valLog.add_error('Missing Sheet',key,f"XLSX {os.path.basename(xlsFile)}",f"Correct XLSX Sheets {list(Submission_Sheets)}")
        fXlsx.close()
        
    return(Submission_Sheets)

#-----------------------------------------------------------------------------
def validate_smiles(smiles_string, *args, **kwargs):
#-----------------------------------------------------------------------------
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    try:
        mol = Chem.MolFromSmiles(smiles_string)
        return mol is not None
    except Exception:
        if valLog:
            valLog.add_error('Wrong SMILES',smiles_string,'Smiles not valid','Correct Smiles')
        return False

#-----------------------------------------------------------------------------
def parse_SampleSubmission(xSheet, valLog=None):
#-----------------------------------------------------------------------------
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
        
        #Smiles
        if _Sample['smiles']:
            try:
                _Sample['smol'] = Chem.MolFromSmiles(_Sample['smiles'])
            except Exception:
                _Sample['smol'] = None
                if valLog:
                    valLog.add_error('Wrong SMILES',_Sample['smiles'],'Smiles not valid','Correct Smiles')
                return False
                
        _Samples.append(_Sample)

    return(_Samples)

#-----------------------------------------------------------------------------
def parse_SampleSubmission(xSheet, valLog=None):
#-----------------------------------------------------------------------------
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
    
    CONTACT_TYPES = {
        0:'PI',
        1:'PC',
        2:'M',        
    }

    _Contacts = []            
    for c in range(3):
        _Contact = {'type':CONTACT_TYPES[c]}
        
        # Col C(2), F(5), I(8)
        for i in CONTACT_FIELDS.keys():
            _Contact[CONTACT_FIELDS[i]] = xSheet[(c*3)+2][i+2]
        
        # Check Email
        try:
            validate_email(_Contact['email'])
        except ValidationError as e:
            if valLog:
                valLog.add_error("Not EMail",_Contact['email'],f" Email of [{_Contact['type']}] not valid","Correct Contact Sheet")   

        _Contacts.append(_Contact)
    return(_Contacts)

#-----------------------------------------------------------------------------
def parse_ContactSubmission(SubmissionFile, SubType=['Samples','Contacts'], **kwargs):
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