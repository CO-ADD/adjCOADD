"""
General base view class or function used by all applications
"""
import os
# from pathlib import Path
from datetime import datetime
from asgiref.sync import sync_to_async

from django import forms
from django.db import transaction, IntegrityError
from django.shortcuts import HttpResponse, render, redirect, get_object_or_404
from django.http import JsonResponse
# from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views import View
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
# from django.conf import settings


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




