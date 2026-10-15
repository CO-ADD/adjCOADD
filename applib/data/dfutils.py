#
import os
import pandas as pd

def sort_pivtable_bylevel(df,nLevel=0):
    _code = list(set([c[nLevel] for c in df.columns]))
    _order = df.columns.reindex(_code, level=0)
    return(df.reindex(columns=_order[nLevel]))

# Load XLSX Sheets into Dictionary 
# -------------------------------------------------
#-----------------------------------------------------------------------------
def get_Xlxs_Sheets(DirName,XlsxFile, SheetDict={}, SheetList=None, FillNA='-', UpperCase=False, **kwargs):
# --------------------------------------------------------------------------------
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
        
    if os.path.isfile(os.path.join(DirName,XlsxFile)):
        fXlsx = open(os.path.join(DirName,XlsxFile), "rb")
        xls = pd.ExcelFile(fXlsx)

        if SheetList is None:
            SheetList = list(SheetDict)
        for key in SheetList:
            if key in xls.sheet_names:
                SheetDict[key] = pd.read_excel(xls, key)
                valLog.add_info('Reading XLSX Sheet',f"{XlsxFile} [{key}]","")
                if UpperCase:
                    SheetDict[key].columns = [c.upper() for c in SheetDict[key].columns]
                else:
                    SheetDict[key].columns = [c.lower() for c in SheetDict[key].columns]
                if FillNA:
                    SheetDict[key] = SheetDict[key].fillna(FillNA)

            else:
                if valLog:
                    valLog.add_error('Missing Sheet',key,f"XLSX {XlsxFile}",f"Correct XLSX Sheets: {SheetList}")
        fXlsx.close()
    
