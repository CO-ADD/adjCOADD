"""
General base view class or function used by all applications
"""
import os
import pandas as pd
import json
from datetime import datetime
from asgiref.sync import sync_to_async

from django import forms
from django.apps import apps
from django.db import transaction, IntegrityError
from django.shortcuts import HttpResponse, render, redirect, get_object_or_404
from django.http import JsonResponse
from django.views import View
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin


# --utilized in Decoration has_permissions, an Alert on Permissions--
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")

# --Super UserRequire Mixin--
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = '/'

    def test_func(self):
        return self.request.user.has_permission('Admin')

# --create view class--
class SimplecreateView(LoginRequiredMixin, View):
    form_class=None
    template_name=None

    def get(self, request, *args, **kwargs):
        form=self.form_class()
        return render(request, self.template_name, {'form':form})
    def post(self, request, *args, **kwargs):
        form =self.form_class(request.POST)
        if form.is_valid():
            with transaction.atomic():
                instance=form.save(commit=False)
                kwargs={'user': request.user}
                instance.save(**kwargs)
            return redirect(request.META['HTTP_REFERER'])
        else:
            messages.error(request, form.errors)
            return redirect(request.META['HTTP_REFERER'])

## create view with searching input:



# --update view class--
class SimpleupdateView(LoginRequiredMixin, View):
    form_class=None
    template_name=None
    model=None

    def get_object_byurlname(self, slug):
        return get_object_or_404(self.model, urlname=slug)
    
    def get_object(self, pk):
        return get_object_or_404(self.model, pk=pk)

    def get(self, request, *args, **kwargs):
        if 'slug' in kwargs:
            slug=kwargs.get("slug")
            object_=self.get_object_byurlname(slug)
        else: 
            pk=kwargs.get("pk")
            object_=self.get_object(pk)
        form=self.form_class(instance=object_)
        return render(request, self.template_name, {'form':form})

    def post(self, request, *args, **kwargs):
        if 'slug' in kwargs:
            slug=kwargs.get("slug")
            object_=self.get_object_byurlname(slug)
        else: 
            pk=kwargs.get("pk")
            object_=self.get_object(pk)
        form =self.form_class(request.POST, instance=object_)
        if form.is_valid():
            print("form is valid")
            with transaction.atomic():
                object_new=form.save(commit=False)
                kwargs={'user': request.user}
                print("form validaed")
                object_new.save(**kwargs)
            return redirect(request.META['HTTP_REFERER'])
        else:
            print("form is not valid")
            messages.error(request, form.errors)
            return redirect(request.META['HTTP_REFERER'])

# --update view class with htmx put request--
from django.http import QueryDict
class HtmxupdateView(LoginRequiredMixin, View):
    form_class=None
    template_name=None
    template_partial=None
    model=None
    
    def get_object(self, pk):
        return get_object_or_404(self.model, pk=pk)

    def get(self, request, *args, **kwargs):
        pk=kwargs.get("pk")
        object_=self.get_object(pk)
        form=self.form_class(instance=object_)
        context={
        "form":form,
        "object":object_,
    }
        return render(request, self.template_name, context)

    def put(self, request, *args, **kwargs):
        pk=kwargs.get("pk")
        object_=self.get_object(pk)
        qd=QueryDict(request.body).dict()
        form =self.form_class(data=qd, instance=object_)
        context={
        "form":form,
        "object":object_,
    }
        if request.GET.get('_value') == 'cancel':
            return render(request, self.template_partial, context)
        elif form.is_valid():
            with transaction.atomic():
                object_new=form.save(commit=False)
                kwargs={'user': request.user}
                object_new.save(**kwargs)                
            return render(request, self.template_partial, context)
        else:
            print("form is not valid")
            messages.error(request, form.errors)
            return render(request, self.template_partial, context)


# export view
class DataExportBaseView(LoginRequiredMixin, View):
    '''
    Export Data in EXCEL or CSV Format
    '''
    def post(self, request):
        app_name = request.POST.get('app_name')
        model_name = request.POST.get('model_name')

        if not app_name or not model_name:
            messages.error(request, "No application or model name were provided for export.")
            return redirect(request.path)

        try:
            Model = apps.get_model(app_name, model_name)
            print(Model)
        except LookupError:
            # Handle the case where the model does not exist.
            return HttpResponse("Model not found.")

        selected_pks_string = request.POST.get('selected_pks')
        print('selected_pks_string:', selected_pks_string)
        if selected_pks_string == 'SelectAll':
            items = Model.objects.all()
        elif selected_pks_string:
            selected_pks = json.loads(selected_pks_string)
            items = Model.objects.filter(pk__in=selected_pks)
        else:
            messages.error(request, "No items were selected for export.")
            return redirect(request.META['HTTP_REFERER'])

        df = pd.DataFrame.from_records(items.values())

        if 'acreated_at' in df.columns:
            df['acreated_at'] = pd.to_datetime(df['acreated_at']).dt.tz_localize(None).apply(lambda a: a.date())
        if 'aupdated_at' in df.columns:
            df['aupdated_at'] = pd.to_datetime(df['aupdated_at']).dt.tz_localize(None).apply(lambda a: a.date())
        if 'adeleted_at' in df.columns:
            df['adeleted_at'] = pd.to_datetime(df['adeleted_at']).dt.tz_localize(None).apply(lambda a: a.date())

        current_date = datetime.now().strftime('%Y-%m-%d')
        filename = f'{model_name}_{current_date}'
        response = HttpResponse("No valid export option was selected.")
        print(request.POST)
        if "csvdownload" in request.POST:
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="{filename}.csv"'
            df.to_csv(path_or_buf=response, index=False)


        elif "xlsxdownload" in request.POST:
            response = HttpResponse(content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
            response['Content-Disposition'] = f'attachment; filename="{filename}.xlsx"'
            df.to_excel(excel_writer=response, index=False)

        return response