import os
import pandas as pd
#
from applib.data.dfutils import get_Xlxs_Sheets


#-----------------------------------------------------------------------------
def get_SeqPrep_xlsx(xlsFile, Sheets=[], FillNA='-', UpperCase=False, **kwargs):
# --------------------------------------------------------------------------------

    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    SeqPrep_Sheets = {
        'Seq': None,
        }

    if os.path.isfile(xlsFile):
        fXlsx = open(xlsFile, "rb")
        xls = pd.ExcelFile(fXlsx)
        #xls = pd.ExcelFile(xlsFile)

        if Sheets is None:
            sheets = [PlatePrep_Sheets.keys()]
        for key in Sheets:
            if key in PlatePrep_Sheets:
                PlatePrep_Sheets[key] = pd.read_excel(xls, key)
                if UpperCase:
                    PlatePrep_Sheets[key].columns = [c.upper() for c in PlatePrep_Sheets[key].columns]
                else:
                    PlatePrep_Sheets[key].columns = [c.lower() for c in PlatePrep_Sheets[key].columns]
                if FillNA:
                    PlatePrep_Sheets[key] = PlatePrep_Sheets[key].fillna(FillNA)
            else:
                if valLog:
                    valLog.add_error('Missing Sheet',key,f"XLSX {os.path.basename(xlsFile)}",f"Correct XLSX Sheets {list(PlatePrep_Sheets)}")
        fXlsx.close()
        
    return(PlatePrep_Sheets)

