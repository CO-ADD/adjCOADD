"""
Used and import for all application filter views
"""

import django_filters
from django import forms
from django.views.generic import ListView
from django.db.models import Q, CharField, TextField, ForeignKey, IntegerField

# -- create a function for search all fields--
def get_all_fields_q_object(model, search_value, exclude_fields=None, prefix=None):
    q_object = Q()
    exclude_fields = exclude_fields or []

    for field in model._meta.get_fields():
        if field.name in exclude_fields:
            continue
        lookup_field_name = f"{prefix}__{field.name}" if prefix else field.name
        if isinstance(field, (CharField, TextField)):
            q_object |= Q(**{f"{lookup_field_name}__icontains": search_value})
        elif isinstance(field, ForeignKey):
            related_model = field.related_model
            related_q_object = get_all_fields_q_object(related_model, search_value, exclude_fields=exclude_fields, prefix=lookup_field_name)
            q_object |= related_q_object
        elif isinstance(field, IntegerField):
            try:
                int_value = int(search_value)
                q_object |= Q(**{lookup_field_name: int_value})
            except ValueError:
                pass
        # Add more field types as needed...

    return q_object
# q_object = get_all_fields_q_object(MyModel, value, exclude_fields=exclude_fields)

# -- Filterset base Class--
class Filterbase(django_filters.FilterSet):
    Search_all_fields = django_filters.CharFilter(method='filter_all_fields', widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder':'Search in All Fields'}),)
    @property
    def qs(self):
        parent = super().qs
        return parent.filter(astatus__gte=0)
   
    def multichoices_filter(self, queryset, name, value):
        lookup='__'.join([name, 'overlap'])
        return queryset.filter(**{lookup: value})
    
    def filter_all_fields(self, queryset, name, value):
        if value:
            exclude_fields = ['password',]
            q_object = get_all_fields_q_object(self._meta.model, value,exclude_fields=exclude_fields)
            return queryset.filter(q_object)
        return queryset


# --Filter view base class--
class FilteredListView(ListView):
    filterset_class = None
    paginate_by=50
    model_fields=None
    order_by=None
    filter_request=None
    filter_Count=None
  
    def get_queryset(self):
        # Get the queryset however you usually would.  For example:
        queryset = super().get_queryset()
        # Then use the query parameters and the queryset to
        # instantiate a filterset and save it as an attribute
        # on the view instance for later.
        self.filterset = self.filterset_class(self.request.GET, queryset=queryset)
        # Return the filtered queryset
        order=self.get_order_by()
        self.filter_Count=self.filterset.qs.distinct().count()
        if order:
            return self.filterset.qs.distinct().order_by(order)
        return self.filterset.qs.distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        self.context_list=context['object_list']
        filter_record_dict={key: self.request.GET.getlist(key) for key in self.request.GET if self.request.GET.getlist(key)!=[""] and key != 'paginate_by' and key!='page' and key!='csrfmiddlewaretoken'}
        filter_record="Selected: "+str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None
        # Pass the filterset to the template - it provides the form.
        context['filter'] = self.filterset
        context['paginate_by']=self.get_paginate_by(self, **kwargs)
        context['fields']=self.model.get_fields(fields=self.model_fields)
        
        context['model_fields']=self.model.get_modelfields(fields=self.model_fields)
        context['filterset']=filter_record
        context['Count']=self.model.objects.count()
        context['querycount']=self.filter_Count
        
        return context

    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)
        return paginate_by
 
    def get_order_by(self):
       
        order_by=self.request.GET.get("order_by", self.order_by) or None
        model_constants_field=self.model_fields
        acs_decs=""
        if order_by:
            order_field=""
            if order_by[0]=="-":
                acs_decs=order_by[0]
                order_field=order_by[1:]
            else:
                order_field=order_by
                
            if order_field in model_constants_field.values():
                order_by=acs_decs+ list(model_constants_field.keys())[list(model_constants_field.values()).index(order_field)]
            elif order_field == 'ID':
                order_by=acs_decs+'pk'
           
            return order_by
        return order_by
