import pandas as pd
import numpy as np
from django.http import JsonResponse, HttpResponse
from django import forms
from django.apps import apps
from django.shortcuts import HttpResponse, render, redirect
from ddrug.models import VITEK_AST, MIC_COADD

# Return styled pivoted table
def flex_pivottable(request,app_model):
    table_html = None

    model_name = app_model.split("-")[1]
    app_name = app_model.split("-")[0]
    pk_list=request.session.get('cached_queryset')
    try:
        Model = apps.get_model(app_name, model_name)
    except LookupError:
        # Handle the case where the model does not exist.
        return HttpResponse("Model not found.")
    model_fields=Model.get_modelfields

    if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
        
        values_str = request.POST.get("values") or "n_left"
        columns_str = request.POST.get("columns") or "orgbatch_id"
        index_str = request.POST.get("index") or "n_created"
        aggfunc_name = request.POST.get("functions") or np.size
        querydata = Model.objects.filter(pk__in=pk_list) if pk_list else Model.objects.all()
        
        values = "n_left" # values_str or None  # pivottable values
        
        if values:
            try:
                # d = {'col1': "n_left", 'col2': pd.Series([2, 3], index=[2, 3])}
                table =  Model.get_pivottable(querydata=querydata, columns_str=columns_str, index_str=index_str, aggfunc=aggfunc_name, values=values)
                #  pd.DataFrame(list(querydata.values())).reset_index()
                # cols = ['n_created','n_left']
                # grbyCol = ['orgbatch_id']
                # df.reindex(columns=['orgbatch_id','n_created','n_left'])
                # table = df[cols].reset_index(False)
                response = HttpResponse(content_type='text/csv')
                response['Content-Disposition'] = 'attachment; filename=pivottable.csv'
                table_html = table.to_html(classes=["table-bordered",])
                table_csv = table.to_csv()

            except Exception as err:
                error_message = f"error is {err}"
                table_html= f"error is {err}"
            return JsonResponse({"table":table_html}) 
    return render(request, 'utils/pivotedtable.html', {"table":table_html, "app_model":app_model, "model_fields":model_fields})