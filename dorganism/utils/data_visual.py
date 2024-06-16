from ddrug.utils.antibiogram import get_Antibiogram_byOrgID
from apputil.utils.data_style import highlight_val, highlight_RSI

## --dOrganism dataframe table, pivottable and styling table--
# convert to dataframe
def convert_dataframe(pk):
    try:
        df = get_Antibiogram_byOrgID(str(pk))
    except Exception as err:
        df = err
        print(f'durgtbl_error {df}')
    return df

# # styling table
# def data_frame_style(pk, displaycols):
#     try:
#         df = convert_dataframe(pk)
#         df.reset_index(inplace=True)
#         df = df[displaycols]
#         df_entries=len(df)
#         style_table=df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
#         table={'df_entries':df_entries, 'style_table': style_table}
#     except Exception as err:
#         table = {'df_entries':err, 'style_table': err}
#     return table

def data_frame_style(pk, displaycols):
    pivdf = get_Antibiogram_byOrgID(str(pk))
    if pivdf:
        pivdf.reset_index(inplace=True)
        pivdf = pivdf[displaycols]
        df_entries=len(pivdf)
        html_table=pivdf.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
        table={'df_entries':df_entries, 'style_table': html_table}
    else:
        table = {'df_entries':0, 'style_table': None}
    return table            
    
# pivottable
# def pivottable_style(pk):
#     '''
#     - using numpy and pandas to pivot table
#     - styling pivottable
#     '''
#     import numpy as np
#     import pandas as pd
#     try:
#         df = convert_dataframe(pk)
#         pivottable = pd.pivot_table(df, columns='BatchID',index=['Drug Class', 'Drug Name', ], values=['BP Profile', 'MIC'],  aggfunc= lambda x:  " ".join([str(y) for y in x]))
#         pivottable = pivottable.astype(str)
#         # Styling pivottable       
#         style_pivot = pivottable.style.applymap(highlight_RSI)     
#         style_pivot = style_pivot.set_table_attributes('class="table table-bordered fixTableHead"').to_html()
#     except Exception as err:
#         style_pivot=err
#     return style_pivot
