from django.http import JsonResponse, HttpResponse
from ddrug.models import VITEK_AST, MIC_COADD
def custom_pivottable(request, querydata, *args, **kwargs):

    selected_data = request.POST.getlist("selected_data[]") or None
    values_str = request.POST.get("values") or None
    columns_str = request.POST.get("columns") or None
    index_str = request.POST.get("index") or None
    aggfunc_name = request.POST.get("functions")
    print(f'values_str is {values_str}')
        
    if selected_data:
        querydata = model.objects.filter(pk__in=selected_data)
         
        
    values = values_str or None  # pivottable values
        
    if values:
        try:
            table = MIC_COADD.get_pivottable(querydata=querydata, columns_str=columns_str, index_str=index_str, aggfunc=aggfunc_name, values=values)
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = 'attachment; filename=pivottable.csv'
            table_html = table.to_html(classes=["table-bordered",])
            table_csv = table.to_csv()
            print(table_html)
            return JsonResponse({"table_html": table_html, "table_csv": table_csv})
        except Exception as err:
            error_message = str(err)
            print(f"error is {err}")
            return JsonResponse({"table_html": error_message,})
    