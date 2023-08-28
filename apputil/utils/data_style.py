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
def highlight_val(val, threhold_number=80.00):
    '''
    highlight the values background
    '''
    if pd.isnull(val): # handle the null values
        val = "empty"
    elif val.replace('.', '', 1).isdigit():  # also consider float values
        val = float(val)
    else:
        val = "empty"
    # background-color settings 
    if val == 'empty':
        bg_color='white'
    elif val > threhold_number:
        bg_color='#dfba9f' 
    else:
        bg_color='yellow'
    
    return 'background-color: %s' % bg_color

# 
def highlight_val2(val):
    '''
    Value Colors,
    R:
    S:
    L:
    '''
    # formatting value 
    if pd.isnull(val):
        val = "empty"
    else:  # also consider float values
        val = str(val)
    # set colors
    if val == 'empty':
        bg_color='white'
    elif val == 'R':
        bg_color='red' 
    elif val == 'S':
        bg_color='green'
    elif val == 'I':
        bg_color='blue'
    elif val == 'nan':
        bg_color = '#FDE992'
    else:
        bg_color='transparent'
    
    return 'background-color: %s' % bg_color
    


def convert_heatmap(xls_file, XlsSheet=None, upload=False, uploaduser='org_db', lower=True):
    df=pd.read_excel(xls_file)
    df.reset_index(drop=True, inplace=True)  # drop the old index
    df = df.astype(str)
    table=df.style.applymap(highlight_val).hide()
    # table=df.style.background_gradient(cmap='Blues')
    table=table.set_table_attributes('class="table table-bordered fixTableHead"') 
    table=table.to_html()
    return table



class Data_visualView(View):
    pass