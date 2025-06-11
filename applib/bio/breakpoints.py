import re
from openpyxl import load_workbook
import pandas as pd

from django.conf import settings

from ddrug.models import Drug, Breakpoint
from apputil.models import ApplicationUser, Dictionary


# ----------------------------------------------------------------------------------------------------
def imp_Breakpoint_fromDict(iDict,valLog,upload=False):
    """
    Create Breakpoint instance from a {Dict}
        drug_name
        org_name
        org_rank
        notorg_name
        notorg_rank
        medical_application     (uncomplicated UTI only) / iv / oral
        bp_resistant_gt
        bp_sensitive_le
        bp_unit                 mg/L / mm
        bp_type                 MIC/ZONE
        bp_source               EUCAST/CLSI
        bp_source_version       v15.0 2025-01-01 / v11.0 2024-12-02
    """
# ----------------------------------------------------------------------------------------------------
    # Change Key names to Lowercase
    iDict =  {k.lower(): v for k, v in iDict.items()} 

    validStatus = True

    DrugID = Drug.get(iDict['drug_name'])
    if DrugID is None:
        validStatus = False
        valLog.add_log('Error','oraOrgDB',f"{iDict['drug_name']} ",'BP Drug does not Exists','-')


    if 'org_name' in iDict:
        OrgName = iDict['org_name']
        OrgRank = Dictionary.get(Breakpoint.DICTIONARY_FIELDS["org_rank"],iDict['org_rank'])
        if OrgRank is None:
            valLog.add_log('Error','oraOrgDB',iDict['org_rank'],'Tax Rank not correct','-')
            validStatus = False
    else:
        OrgName = None
        OrgRank = None

    if 'notorg_name' in iDict:
        NotOrgName = iDict['notorg_name']
        NotOrgRank = Dictionary.get(Breakpoint.DICTIONARY_FIELDS["notorg_rank"],iDict['notorg_rank'])
        if NotOrgRank is None:
            valLog.add_log('Error','oraOrgDB',iDict['notorg_rank'],'(Not) Tax Rank not correct','-')
            validStatus = False
    else:
        NotOrgName = None
        NotOrgRank = None

    djBP = Breakpoint.get(DrugID, OrgName, OrgRank, NotOrgName, NotOrgRank,
                        iDict['medical_application'], iDict['bp_type'], iDict['bp_source'])
    if djBP is None:
        djBP = Breakpoint()
        djBP.drug_id = DrugID
        djBP.org_name = OrgName
        djBP.org_rank = OrgRank
        djBP.notorg_name = NotOrgName
        djBP.notorg_rank = NotOrgRank

        valLog.add_log('Info',"",f"{iDict['drug_name']} {OrgRank} {OrgName} {NotOrgRank} {NotOrgName}",'New BP','-')

    djBP.bp_type = Dictionary.get(djBP.DICTIONARY_FIELDS["bp_type"],iDict['bp_type'])
    if djBP.bp_type is None:
        valLog.add_log('Error','oraOrgDB',iDict['bp_type'],'BP Type not correct','-')
        validStatus = False

    djBP.med_application = iDict['medical_application']
    djBP.bp_res_gt = iDict['bp_resistant_gt']
    djBP.bp_sens_le = iDict['bp_sensitive_le']
    djBP.bp_unit = iDict['bp_unit']
    djBP.bp_comb = iDict['combination_type']
    djBP.bp_source = iDict['bp_source']
    djBP.bp_source_version = iDict['bp_source_version']

    djBP.set_defaults_model()
    validStatus = True
    validDict = djBP.validate_fields()
    if validDict:
        validStatus = False
        for k in validDict:
            valLog.add_log('Warning','',k,validDict[k],'-')
            #print(f"Warning : {k} {validDict[k]}")
    djBP.VALID_STATUS = validStatus

    return(djBP)


#----------------------------------------------------------------------------
class EUCAST():
#----------------------------------------------------------------------------
    BP_COLUMN_NAME  = "BP_INDEX"

    EUCAST_SHEETS = {
        'Enterobacterales':'Enterobacterales',
        'Pseudomonas':'Pseudomonas',
        'S.maltophilia':'Stenotrophomonas maltophilia',
        'Acinetobacter':'Acinetobacter',
        'Staphylococcus':'Staphylococcus',
        'Enterococcus':'Enterococcus',
        'Streptococcus A,B,C,G':'Streptococcus',
        'S.pneumoniae':'Streptococcus pneumoniae',
        'H.influenzae':'Haemophilus influenzae',
        'M.catarrhalis':'Moraxella catarrhali',
        'N.gonorrhoeae':'Neisseria gonorrhoeae',
        'N.meningitidis':'Neisseria meningitidis',
        'H.pylori':'Helicobacter pylori',
        'L.monocytogenes':'Listeria monocytogenes',
        'Pasteurella':'Pasteurella',
        'Corynebacterium':'Corynebacterium',
        'M.tuberculosis':'Mycobacterium tuberculosis',
        'L.pneumophila':'Ligonella pneumophila',
        'B.cepacia':'Burkholderia cepacia',
        'B.pseudomallei':'Burkholderia pseudomallei',
        'B.melitensis ':'Brucella melitensis ',
        'B.anthracis':'Bacillus anthracis',
        'Bacillus':'Bacillus',
        'A.xylosoxidans':'Achromobacter xylosoxidans',
        'Vibrio':'Vibrio',
        }
    
    def __init__(self,XlsName, Version='v15.0 2025-01-01', BPType='MIC', BPClass='Bacteria'):
        self.xlsx_name = XlsName
        self.bp_type = BPType
        self.bp_class = BPClass
        self.bp_source = 'EUCAST'
        self.bp_source_version = Version
        
        self.data = []
        self.df = None

    #----------------------------------------------------------------------------        
    def process_sheets(self, sheet_dict = EUCAST_SHEETS):
    #----------------------------------------------------------------------------        
        self.df = pd.DataFrame()
        for s in sheet_dict:
            print(f" [EUCAST-Sheet] {s} - OrgName: {sheet_dict[s]}")
            _df = self.process_single_sheet(s)
            _df['ORG_NAME'] = sheet_dict[s]
            _df['BP_SOURCE'] = self.bp_source
            _df['BP_SOURCE_VERSION'] = self.bp_source_version
            
            self.df = pd.concat([self.df,_df])


    #----------------------------------------------------------------------------        
    def process_single_sheet(self, sheetname):
    #----------------------------------------------------------------------------        
        _df = self._load_raw_table(sheetname)
        _df = self._select_bp_data(_df)
        return _df

    #----------------------------------------------------------------------------        
    def _load_raw_table(self,sheet_name):
    #----------------------------------------------------------------------------        

        # -- Load workbook -------------------------------------------
        wb = load_workbook(self.xlsx_name, rich_text=True, data_only=True)

        # Check if the sheet_name exists
        if sheet_name not in wb.sheetnames:
            raise ValueError(f"Sheet '{sheet_name}' not found in the workbook.")

        # Process only the matching sheet
        sheet = wb[sheet_name]
        data = []

        # -- Get Raw datatable -----------------------------------------
        for row in sheet.iter_rows():
            #row_dict = {'Sheet':sheet_name}
            row_data = [sheet_name]
            for cell in row:
                if isinstance(cell.value, str):
                    cleaned_value = cell.value
                    cleaned_value = re.sub(r"<.*?>", "", cleaned_value)  # Remove HTML-like tags
                    cleaned_value = re.sub(r"[\']", "", cleaned_value)  # Remove '
                    #row_dict['drug'] = cleaned_value
                    row_data.append(cleaned_value)
                    
                elif isinstance(cell.value, list):
                    cleaned_value = cell.value[0]
                    if isinstance(cleaned_value, str) and cleaned_value == "":
                        cleaned_value = cell.value[1]
                    row_data.append(str(cleaned_value))
                else:
                    row_data.append(str(cell.value))  # Append non-string values as-is
            data.append(row_data)  # Append the processed row
            
        # Return the processed DataFrame for the matching sheet
        return(pd.DataFrame(data))
    
    #----------------------------------------------------------------------------        
    @staticmethod
    def apply_drugname(s):
    #----------------------------------------------------------------------------        
        s['MEDICAL'] = ""
        
        _n = s["DRUG_NAME"].strip()  
        if _n:
            # Find iv/oral
            for _m in ["iv","oral"]:
                if _m in _n:
                    _n = _n.replace(_m,"").strip()
                    s['MEDICAL'] = f"{s['MEDICAL']} {_m}"

            # Find ()
            _r = re.search(r'\((.*?)\)', _n)
            if _r:
                s['MEDICAL'] = f"{s['MEDICAL']} - {_r.group(1)}"
                _n = _n.replace(f"({_r.group(1)})","").strip()
            
            # Repace combination.sep
            _s = _n.split('-')
            if len(_s) > 1:
                _n = '|'.join([_s[0].capitalize(),_s[1].capitalize()])
        
            # Final DrugName
            s["DRUG_NAME"] = _n
            
        # Reset BP_INDEX
        if s['BP_MIC_R_GT'] == '-' and s['BP_ZONE_R_GT'] == '-':
            s['BP_INDEX'] = 0
        
        if s['BP_MIC_R_GT'] == 'None' and s['BP_ZONE_R_GT'] == 'None':
            s['BP_INDEX'] = 0

        return(s)
    
    #----------------------------------------------------------------------------        
    def _select_bp_data(self,data_df):
    #----------------------------------------------------------------------------        
    
        _n_col = len(data_df.columns)
        print(data_df.head())
        if _n_col < 10:
            # For Sheest with only MIC (no ZONE)
            data_df.columns = ['SHEET','DRUG_NAME',
                            'BP_MIC_R_GT','BP_MIC_S_LE','BP_MIC_ATU',
                            'NOTES',
                            #self.BP_COLUMN_NAME
                            ] + list(data_df.columns[6:])
            
            for c in ['BP_ZONE_CONC','BP_ZONE_R_GT','BP_ZONE_S_LE','BP_ZONE_ATU']:
                data_df[c] = '-'

        else:    
            data_df.columns = ['SHEET','DRUG_NAME',
                            'BP_MIC_R_GT','BP_MIC_S_LE','BP_MIC_ATU',
                            'BP_ZONE_CONC','BP_ZONE_R_GT','BP_ZONE_S_LE','BP_ZONE_ATU',
                            'NOTES',
                            #self.BP_COLUMN_NAME
                            ] + list(data_df.columns[10:])
        #-- Extract MIC/ZONE data -----------------------------------------
        # Preprocessing
        data_df = data_df[data_df['DRUG_NAME'].notnull()]
        data_df = data_df[data_df['DRUG_NAME'] != 'None']
        
        # 2 - Column with "MIC breakpoint"
        data_df.loc[:, 'BP_MIC_R_GT'] = data_df['BP_MIC_R_GT'].astype(str)
        data_df[self.BP_COLUMN_NAME] = data_df['BP_MIC_R_GT'].str.find("MIC breakpoint")

        inside_mic = False

        # Search for MIC 
        for idx in data_df.index:
            mic_entry = data_df.at[idx, self.BP_COLUMN_NAME]
            if mic_entry == 0:  # Found the 'MIC breakpoint'
                inside_mic = True
                continue
            if inside_mic:
                if mic_entry == -1:  # No 'MIC breakpoint'
                    data_df.at[idx, self.BP_COLUMN_NAME] = 1.0
                if pd.isna(
                    data_df.at[idx, self.BP_COLUMN_NAME]
                ):  # Reset when a row is NaN (indicating end of MIC section)
                    inside_mic = False

        # Filter for MIC data
        data_df = data_df[data_df[self.BP_COLUMN_NAME] == 1.0]
        
        # Remove '()'
        for col in ['BP_MIC_R_GT','BP_MIC_S_LE','BP_ZONE_CONC','BP_ZONE_R_GT','BP_ZONE_S_LE']:
            data_df[col] = data_df[col].map(lambda x: re.sub(r"[()]", "", x) if isinstance(x, str) else x)            
        
        # Fix DrugName
        data_df = data_df.dropna(subset=["DRUG_NAME"])
        data_df["DRUG_NAME"] = data_df["DRUG_NAME"].map(lambda x: re.sub(r",", "", x) if isinstance(x, str) else x)
        data_df = data_df.apply(self.apply_drugname,axis=1)
        data_df = data_df[data_df['BP_INDEX'] > 0]
        
        return data_df

