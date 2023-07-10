import ddrug.utils.tables as drugtbl
from apputil.utils.data_style import highlight_val, highlight_val2

## --dOrganism dataframe table, pivottable and styling table--
# convert to dataframe
def convert_dataframe(pk):
    print(pk)
    try:
        df = drugtbl.get_Antibiogram_byOrgID(pk)
        print(f'dataframe: {df}')
    except Exception as err:
        df = err
        print(f'error {df}')
   
    return df

# styling table
def data_frame_style(pk, displaycols):
    try:
        df = convert_dataframe(pk)
        df.reset_index(inplace=True)
        df = df[displaycols]
        df_entries=len(df)
        style_table=df.to_html(classes=["dataframe", "table", "table-bordered", "fixTableHead"], index=False)
        table={'df_entries':df_entries, 'style_table': style_table}
    except Exception as err:
        table = {'df_entries':err, 'style_table': err}
    return table

# pivottable
def pivottable_style(pk):
    '''
    - using numpy and pandas to pivot table
    - styling pivottable
    '''
    import numpy as np
    import pandas as pd
    try:
        df = convert_dataframe(pk)
        pivottable = pd.pivot_table(df, columns='BatchID',index=['Drug Class', 'Drug Name', ], values=['BP Profile', 'MIC'],  aggfunc= lambda x:  " ".join([str(y) for y in x]))
        pivottable = pivottable.astype(str)
        # Styling pivottable       
        style_pivot = pivottable.style.applymap(highlight_val2)     
        style_pivot = style_pivot.set_table_attributes('class="table table-bordered fixTableHead"').to_html()
    except Exception as err:
        style_pivot=err
    return style_pivot
