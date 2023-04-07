from rest_framework import generics
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework.pagination import PageNumberPagination
# from .models import MyModel
# from .serializers import MyModelSerializer
# from .filtersets import MyModelFilterSet

class CustomPagination(PageNumberPagination):
    page_size = 20
    page_size_query_param = 'page_size'
    max_page_size = 1000


class API_FilteredListView(generics.ListAPIView):
    queryset = None
    serializer_class = None
    filterset_class = None
    filter_backends = [DjangoFilterBackend]
    # paginate_by = 20
    model_fields = None
    order_by = None
    filter_request=None
    filter_Count=None
    pagination_class = CustomPagination

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

    def get_serializer_context(self):
        context = super().get_serializer_context()
        context['filter'] = self.filterset
        return context

    # def get_paginate_by(self, queryset):
    #     paginate_by = self.request.GET.get("paginate_by", self.paginate_by)
    #     return paginate_by

    def get_order_by(self):
        order_by = self.request.GET.get("order_by", self.order_by) or None
        model_constants_field = self.model_fields
        acs_decs = ""
        if order_by:
            order_field = ""
            if order_by[0] == "-":
                acs_decs = order_by[0]
                order_field = order_by[1:]
            else:
                order_field = order_by

            if order_field in model_constants_field.values():
                order_by = acs_decs + list(model_constants_field.keys())[list(model_constants_field.values()).index(order_field)]
            elif order_field == 'ID':
                order_by = acs_decs + 'pk'

        return order_by