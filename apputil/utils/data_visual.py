"""
Data Visulization functions:
- Heat Map
"""

import pandas as pd
from django.views import View
from django.shortcuts import render

def highlight_val(val):
    '''
    highlight the value
    '''
    if pd.isnull(val):
        val = "empty"
    else:
        val =str(val)
    bg_color= 'white' if val == 'empty' else 'var(--bs-inform-color)'
    return 'background-color: %s' % bg_color


def convert_heatmap(xls_file, XlsSheet=None, upload=False, uploaduser='org_db', lower=True):
    print(xls_file)
    print(XlsSheet)
    df=pd.read_excel(xls_file)
    df.reset_index(inplace=True)
    table=df
    table=table.style.applymap(highlight_val).set_table_attributes('class="table table-bordered fixTableHead"') 
    table=df.to_html()
    return table

class Data_visualView(View):
    pass