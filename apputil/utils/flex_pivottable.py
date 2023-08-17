import pandas as pd
import numpy as np
from django.apps import apps
from django.shortcuts import HttpResponse, render
from apputil.utils.data_style import highlight_val, highlight_val2
from ddrug.utils.bio_data import agg_Lst, agg_DR, agg_Inhib

#    #-------------------------------------------------------------------------------------------------
# data-visulization 
def get_pivottable(querydata=None, aggfunc_table=None, columns_table=None, index_table=None, values=None):
    pivot_table={}
    aggfunc={'String_concat':lambda x:  " ".join([str(y) for y in x]), 'Lst': agg_Lst, 'DR': agg_DR, 'Inhib':agg_Inhib}
    data=querydata #
    df=pd.DataFrame(data)
    df.reset_index
    df.fillna(0)
   
    columns=columns_table.split(",") 
    index=index_table.split(",")
  
    try:
        table=pd.pivot_table(df, values=values, index=index,
                       columns=columns, aggfunc=aggfunc[aggfunc_table], fill_value=0)#np_aggfunc[aggfunc], fill_value='0')

    except Exception as err:
        print(f'err is {err}')
        table=err
    # pivot_table={columns:columns, index: index, table: table}
    pivot_table["columns"]=columns
    pivot_table["index"]=index
    pivot_table["table"]=table
    return pivot_table
# Return styled pivoted table
def flex_pivottable(request, app_model):
    table_html = None
    table = None
    model_name = app_model.split("-")[1]
    app_name = app_model.split("-")[0]
    values_str = None
    select_hrfields = None
    select_vtfields = None
    select_fields = None
    select_value = None

    try:
        Model = apps.get_model(app_name, model_name)
    except LookupError:
        # Handle the case where the model does not exist.
        return HttpResponse("Model not found.")
    model_fields=[f.replace(".", "__") for f in list(Model.HEADER_FIELDS.keys())] #list(Model.HEADER_FIELDS.keys())
    pk_list = request.session.get(f'{Model}_cached_queryset')
     # get queryset
    queryset = Model.objects.filter(pk__in=pk_list) if pk_list else Model.objects.all()
    # flatten the queryset to put it in dataframe
    querylist=queryset.values_list(*model_fields)
    data={}
    for num_fields in range(len(model_fields)):
        arrary=[querylist[i][num_fields] for i in range(len(querylist))]
        data[model_fields[num_fields]]=arrary

    if "download" in request.POST:
        values_table = request.session.get(f"{request.user}_pivoteddata")[0]
        columns_table = request.session.get(f"{request.user}_pivoteddata")[1]
        index_table = request.session.get(f"{request.user}_pivoteddata")[2]
        aggfunc_table = request.session.get(f"{request.user}_pivoteddata")[3]

        table = get_pivottable(querydata=data, aggfunc_table = aggfunc_table, columns_table = columns_table, index_table= index_table, values = values_table)['table']
        response = HttpResponse(content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = f'attachment; filename="pivot.xlsx"'
        table_excel = table.to_excel(excel_writer=response, index=True)
        return response

    if request.method == "POST": # 
        
        values_str = request.POST.get("values") or None  # or "n_left"
        columns_str =request.POST.get("column") or None # or "orgbatch_id"
        index_str = request.POST.get("index") or None 
        aggfunc_name = request.POST.get("data_function_type") or None
        request.session[f"{request.user}_pivoteddata"]= [values_str, columns_str, index_str, aggfunc_name]
       
        values = values_str or None  # pivottable values
        
        if values:
            select_value = values
            try:    
                result = get_pivottable(querydata=data, aggfunc_table=aggfunc_name, columns_table=index_str, index_table=columns_str, values=values_str)
                table = result["table"]
    
                select_hrfields = result["columns"]
                select_vtfields = result["index"]
                select_fields = select_hrfields + select_vtfields
                try:
                    # table_html = table.astype(str)
                    # print("html")
                    # table_html = table.style.applymap(hightlight_val2)
                    # print("set attr")
                    # table_html= table_html.set_table_attributes('class="table table-bordered fixTableHead"').to_html()
                    table_html = table.head(n=10).reindex().to_html(classes=["table", "fixTableHead", "overflow-auto", "table-hover"])
                except Exception as err:
                    table_html = f"<span class='text-danger'>something wrong with {table}</span>"
 
            except Exception as err:
                error_message = f"error is {err}"
                table_html= f"error is {err}"
        
        
       

    return render(request, 'utils/pivotedtable.html', {"table":table_html, "select_value": select_value, "select_fields":select_fields, "select_hrfields":select_hrfields, "select_vtfields":select_vtfields, "app_model":app_model, "model_fields":model_fields, "query_list": pk_list})
