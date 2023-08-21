import pandas as pd
import numpy as np
from django.apps import apps
from django.shortcuts import HttpResponse, render
from apputil.utils.data_style import highlight_val, highlight_val2
from ddrug.utils.bio_data import agg_Lst, agg_DR, agg_Inhib

#    #-------------------------------------------------------------------------------------------------
# get_pivottable: - used to convert queryset into pivottable 
#                 - return pivotedtable with user selected columns, index and values
def get_pivottable(querydata=None, aggfunc_table=None, columns_table=None, index_table=None, values=None, fields_dict = None):
    pivot_table = {}
    aggfunc = {'String_concat':lambda x:  " ".join([str(y) for y in x]), 'Lst': agg_Lst, 'DR': agg_DR, 'Inhib':agg_Inhib}
    # selected data, columns, index, value from UI
    data_ui = querydata 
    columns_ui = columns_table.split(",")
    columns_pivot = [list(fields_dict[col].keys())[0] if isinstance(fields_dict[col], dict) else fields_dict[col] for col in columns_ui]
    index_ui=index_table.split(",")
    index_pivot = [list(fields_dict[i].keys())[0] if isinstance(fields_dict[i], dict) else fields_dict[i] for i in index_ui]
    # generate dataframe with selected data,
    # rename columns' name to verbose name via fields_dict
    df = pd.DataFrame(data_ui)  
    table_columns = {col: list(fields_dict[col].keys())[0] if isinstance(fields_dict[col], dict) else fields_dict[col]  for col in df.columns }
    df = df.rename(columns = table_columns)
    table_values = list(fields_dict[values].keys())[0] if isinstance(fields_dict[values], dict) else fields_dict[values]
    try:
        table=pd.pivot_table(df, values=[table_values,], index=index_pivot,
                       columns=columns_pivot, aggfunc=aggfunc[aggfunc_table],)
    except Exception as err:
        print(f'err is {err}')
        table = err
    pivot_table={'columns':columns_ui, 'index': index_ui, 'table': table, 'table_values': table_values}

    return pivot_table

# flex_pivottable used to: - process and pass user selected data send to func get_pivottable,
#                          - convert modelfield name to verbose name 
#                          - recieve pivotted table and columns, index from func,
#                          - to style table and send to UI
def flex_pivottable(request, app_model):
    # table = None
    # get model from UI
    model_name = app_model.split("-")[1]
    app_name = app_model.split("-")[0]
    table_html = None # render styled table to UI
    fields_dict = {} # fields_dict contains key(modelfields): modelfields' verbosename
    values_str = None # model field
    select_hrfields = None # model field's verbosename
    select_vtfields = None
    select_fields = None
    select_value = None

    try:
        Model = apps.get_model(app_name, model_name)
    except LookupError:
        # Handle the case where the model does not exist.
        return HttpResponse("Model not found.")
    model_fields = [f.replace(".", "__") for f in list(Model.HEADER_FIELDS.keys())] 
    verbose_fields = Model.get_fields()
    for i in range(len(model_fields)):
        fields_dict = {model_fields[i]: verbose_fields[i] for i in range(len(verbose_fields)) }
    pk_list = request.session.get(f'{Model}_cached_queryset')
    # get queryset
    queryset = Model.objects.filter(pk__in=pk_list) if pk_list else Model.objects.all()
    # flatten the queryset to put it in dataframe
    querylist=queryset.values_list(*model_fields)
    data={}
    for num_fields in range(len(model_fields)):
        arrary=[querylist[i][num_fields] for i in range(len(querylist))]
        data[model_fields[num_fields]]=arrary
    #
    #  When Click download button 
    if "download" in request.POST:
        values_table_dl = request.session.get(f"{request.user}_pivoteddata")[0]
        columns_table_dl = request.session.get(f"{request.user}_pivoteddata")[1]
        index_table_dl = request.session.get(f"{request.user}_pivoteddata")[2]
        aggfunc_table_dl = request.session.get(f"{request.user}_pivoteddata")[3]
        result = get_pivottable(querydata=data, aggfunc_table=aggfunc_table_dl, columns_table = index_table_dl, index_table = columns_table_dl, values=values_table_dl, fields_dict=fields_dict)
        table = result["table"]
        response = HttpResponse(content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = f'attachment; filename="pivot.xlsx"'
        table_excel = table.to_excel(excel_writer=response, index=True)
        return response
    # 
    #  When click go arrow button (Generate Pivot table)
    if request.method == "POST": #       
        values_str = request.POST.get("values") or None  # or "n_left"
        columns_str =request.POST.get("column") or None # or "orgbatch_id"
        index_str = request.POST.get("index") or None 
        aggfunc_name = request.POST.get("data_function_type") or None
        request.session[f"{request.user}_pivoteddata"]= [values_str, columns_str, index_str, aggfunc_name]
        values = values_str or None  # pivottable values
        
        if values:
            try:    
                result = get_pivottable(querydata=data, aggfunc_table=aggfunc_name, columns_table=index_str, index_table=columns_str, values=values_str, fields_dict=fields_dict)
                table = result["table"]
                table = table.head(n=10)  \
                             .astype(str) \
                             .style       \
                             .applymap(highlight_val2) \
                             .set_table_attributes('class="table table-bordered fixTableHead"')

                select_value = result["table_values"]
                select_hrfields ={result["columns"][i]: fields_dict[result["columns"][i]] for i in range(len(result["columns"]))} #result["columns"]
                select_vtfields ={result["index"][i]: fields_dict[result["index"][i]] for i in range(len(result["index"]))} # result["index"]
                select_fields = list(select_hrfields.keys()) + list(select_vtfields.keys())
                try:
                    table_html = table.to_html()
                except Exception as err:
                    table_html = f"<span class='text-danger'>something wrong with {table}</span>"
 
            except Exception as err:
                error_message = f"error is {err}"
                table_html= f"error is {err}"
        
    return render(request, 'utils/pivotedtable.html', {"table":table_html, "select_value": select_value, "select_fields":select_fields, "select_hrfields":select_hrfields, "select_vtfields":select_vtfields, "app_model":app_model, "model_fields":fields_dict, "query_list": pk_list})
