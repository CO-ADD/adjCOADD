"""
Data Visulization functions:
- Heat Map
"""

import pandas as pd
import numpy as np
from django.views import View
from django.shortcuts import render
from .files_upload import file_location

## Different stylings for dataframe table, pivottable, etc...
# Highlighted Table styles
def highlight_val(val):
    '''
    highlight the values background
    '''
    if pd.isnull(val):
        val = "empty"
    elif val.replace('.', '', 1).isdigit():  # also consider float values
        val = float(val)
    else:
        val = "empty"
    if val == 'empty':
        bg_color='white'
    elif val > 80.00:
        bg_color='#dfba9f' 
    else:
        bg_color='transparent'
    
    return 'background-color: %s' % bg_color

# 
def hightlight_val2(val):
    pass
    


def convert_heatmap(xls_file, XlsSheet=None, upload=False, uploaduser='org_db', lower=True):
    df=pd.read_excel(xls_file)
    df.reset_index(drop=True, inplace=True)  # drop the old index
    df = df.astype(str)
    table=df.style.applymap(highlight_val).hide_index()
    # table=df.style.background_gradient(cmap='Blues')
    table=table.set_table_attributes('class="table table-bordered fixTableHead"') 
    table=table.to_html()
    print(table)
    return table



class Data_visualView(View):
    pass