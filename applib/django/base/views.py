"""
General base view class or function used by all applications
"""
import os
import pandas as pd
import json
from datetime import datetime

from django.apps import apps
from django.core.files.storage import default_storage
from django.core.exceptions import ValidationError
from django.db import transaction, IntegrityError
from django.shortcuts import HttpResponse, render, redirect, get_object_or_404
from django.http import JsonResponse, QueryDict
from django.views import View
from django.views.generic import ListView
from django.views.generic.edit import FormView
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin

from ddrug.utils.antibiogram import get_Antibiogram_byOrgID,piv_Antibiogram_byOrgID
from apputil.models import ApplicationLog
from apputil.utils.files_upload import SuperUserRequiredMixin


# -----------------------------------------------------------------
# --utilized in Decoration has_permissions, an Alert on Permissions--
# -----------------------------------------------------------------

def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")

# # -----------------------------------------------------------------
# # --Super UserRequire Mixin--
# # -----------------------------------------------------------------
# class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
#     login_url = '/'

#     def test_func(self):
#         return self.request.user.has_permission('Admin')

#     def handle_no_permission(self):
#         return HttpResponse( 'Only users with ADMIN permission have access to this view')

# -----------------------------------------------------------------
# --Write UserRequire Mixin--
# -----------------------------------------------------------------
class WriteUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):
    login_url = '/'

    def test_func(self):
        return self.request.user.has_permission('Write')
    
    def handle_no_permission(self):
        return HttpResponse( 'Only users with WRITE permission have access to this view')

# -----------------------------------------------------------------
# -- Create View class--
# -----------------------------------------------------------------
class Base_CreateView(LoginRequiredMixin, View):
    form_class = None
    template_name = None
    transaction_use = 'default'

    # -----------------------------------------------------------
    def get(self, request, *args, **kwargs):
        form=self.form_class()
        return render(request, self.template_name, {'form':form})
    
    # -----------------------------------------------------------
    def post(self, request, *args, **kwargs):
        
        form =self.form_class(request.POST, request.FILES)
        if form.is_valid():
            with transaction.atomic(using=self.transaction_use):
                instance=form.save(commit=False)
                kwargs={'user': request.user}
                instance.save(**kwargs)
                ## python logging levels: 10-'DEBUG', 40-'ERROR', 50-'CRITICAL', 30-'WARNING', 20-'INFO', 0-'Notset'
                 #LogCode, LogProc,LogType,LogUser,LogObject,LogDesc,LogStatus
                ApplicationLog.add('Create',str(instance.pk),'Info',request.user,str(instance.pk)[:10],'Create a new entry','Completed')
            return redirect(request.META['HTTP_REFERER'])
        else:
            messages.error(request, form.errors)
            return redirect(request.META['HTTP_REFERER'])

# -----------------------------------------------------------------
# -- Update View class --
# -----------------------------------------------------------------
class Base_UpdateView(LoginRequiredMixin, View):
    form_class = None
    template_name = None
    model = None
    transaction_use = 'default'

    # -----------------------------------------------------------
    def get_object_byurlname(self, slug):
        return get_object_or_404(self.model, urlname=slug)
    
    # -----------------------------------------------------------
    def get_object(self, pk):
        return get_object_or_404(self.model, pk=pk)

    # -----------------------------------------------------------
    def get(self, request, *args, **kwargs):
        if 'slug' in kwargs:
            slug=kwargs.get("slug")
            object_=self.get_object_byurlname(slug)
        else: 
            pk=kwargs.get("pk")
            object_=self.get_object(pk)
            print("test")
        form=self.form_class(instance=object_)
        return render(request, self.template_name, {'form':form})

    # -----------------------------------------------------------
    def post(self, request, *args, **kwargs):
        if 'slug' in kwargs:
            slug=kwargs.get("slug")
            object_=self.get_object_byurlname(slug)
        else: 
            pk=kwargs.get("pk")
            object_=self.get_object(pk)
        form =self.form_class(request.POST, instance=object_)
        if form.is_valid():
            with transaction.atomic(using=self.transaction_use):
                object_new=form.save(commit=False)
                kwargs={'user': request.user}
                object_new.save(**kwargs)
                ApplicationLog.add('Update',str(object_new.pk),'Info',request.user,str(object_new.pk),'Update an entry','Completed')
            return redirect(request.META['HTTP_REFERER'])
        else:
            messages.error(request, form.errors)
            return redirect(request.META['HTTP_REFERER'])

# -----------------------------------------------------------------
# -- Remove View class --
# -----------------------------------------------------------------
class Base_RemoveView(SuperUserRequiredMixin, Base_UpdateView):
    model=None
    transaction_use = 'default'
    
    # -----------------------------------------------------------
    def post(self, request, *args, **kwargs):
        if 'slug' in kwargs:
            slug=kwargs.get("slug")
            object_=self.get_object_byurlname(slug)
        else: 
            pk=kwargs.get("pk")
            object_=self.get_object(pk)
        with transaction.atomic(using=self.transaction_use):
            kwargs={'user': request.user}
            try:
                object_.remove(**kwargs)
                ApplicationLog.add('Removed','log_proc','Warning',request.user, str(object_.pk), 'switch entry_astatus -9','Completed')            
            except Exception as err:
                messages.error(request, err)

        return redirect(request.META['HTTP_REFERER'])
    
# -----------------------------------------------------------------
# --update view class with htmx put request--
# -----------------------------------------------------------------
class Htmx_UpdateView(LoginRequiredMixin, View):
    form_class = None
    template_name = None
    template_htmx = None
    model = None
    transaction_use = 'default'
    
    def get_object(self, pk):
        return get_object_or_404(self.model, pk=pk)

    def get(self, request, *args, **kwargs):
        pk=kwargs.get("pk")
        object_=self.get_object(pk)
        form=self.form_class(instance=object_)
        context={"form":form,
                 "object":object_,
                }
        return render(request, self.template_name, context)

    def put(self, request, *args, **kwargs):
        pk=kwargs.get("pk")
        object_=self.get_object(pk)
        qd=QueryDict(request.body).dict()
        form =self.form_class(data=qd, instance=object_)
        context={"form":form,
                "object":object_,
                }        
        if request.GET.get('_value') == 'cancel':
            return render(request, self.template_htmx, context)
        
        elif form.is_valid():
            with transaction.atomic(using=self.transaction_use):
                object_new=form.save(commit=False)
                kwargs={'user': request.user}
                object_new.save(**kwargs)
                ApplicationLog.add('Update',str(object_new.pk),'Info', request.user, str(object_new.pk),f'Update an {object_new._meta.model}','Completed')              
                #print(f'Update an {object_new._meta.model} : {str(object_new.pk)}')
            return render(request, self.template_htmx, context)
        
        else:
            # raise ValidationError
            context["form_errors"] = form.errors
            # messages.error(request, form.errors)
            return render(request, self.template_htmx, context)


# -----------------------------------------------------------------
# -- View class for Filtered List
# -----------------------------------------------------------------
class Filtered_ListView(ListView):
    
    filterset_class = None #each filterset class based on class BaseStatus_Filter
    paginate_by = 50
    model_fields = None
    order_by = None
    filter_Count = None
    app_name = None
    model_name = None

    #--------------------------------------------------------------------------
    @staticmethod
    def find_item_index(lst, item):
        for i, element in enumerate(lst):
            if isinstance(element, dict):
                if item in element.keys():
                    return i
            elif element == item:
                return i
        return -1
  
    #--------------------------------------------------------------------------
    def get_filter_record(self):
        EXCLUDE_KEYS = ['paginate_by','page', 'csrfmiddlewaretoken', 'reset', "pivot", "applysingle", "applymulti"]
    
        # filter_record_dict = {}
        # _Excluded_Keys = ['paginate_by','page', 'csrfmiddlewaretoken', 'reset', "pivot", "applysingle", "applymulti"]
        # for key in self.request.GET:
        #     if key not in EXCLUDE_KEYS:
        #         if self.request.GET.getlist(key)!=[""] :
        #             filter_record_dict[key] = self.request.GET.getlist(key)
    
        return({key: self.request.GET.getlist(key) for key in self.request.GET if self.request.GET.getlist(key)!=[""] and key not in EXCLUDE_KEYS})
  
    #--------------------------------------------------------------------------
    def get_queryset(self):
 
        # Get the queryset however you usually would.  For example:
        queryset = super().get_queryset()
        kwargs={'deep': False}      

        # Check if the reset request is submitted
        # Remove the stored queryset from the session
        if self.request.GET.get('reset')=='True':
            if 'cached_queryset' in self.request.session:
                del self.request.session[f'{self.model}_cached_queryset'] 
                
        # Instantiate the filterset with either the stored queryset from the session or the default queryset
        # ---- Switch off cache queryset
        # if self.request.session.get('cached_queryset'):
        #     stored_queryset_pks = self.request.session['cached_queryset']
        #     stored_queryset = queryset.filter(pk__in=stored_queryset_pks)
        # ----
        
        filter_record_dict = self.get_filter_record()
                        
        if 'applymulti' in self.request.GET:
            kwargs={'deep': True}
            self.filterset = self.filterset_class(self.request.GET,  queryset = queryset, filterset_dict= filter_record_dict, **kwargs)
        else:
            self.filterset = self.filterset_class(self.request.GET, queryset = queryset, filterset_dict= filter_record_dict, **kwargs)
            
        # Cache the filtered queryset in the session
        filtered_queryset_pks = self.filterset.qs.distinct().values_list('pk', flat = True)
        self.request.session[f'{self.model}_cached_queryset'] = list(filtered_queryset_pks) if filtered_queryset_pks else None  

        # Then use the query parameters and the queryset to
        # instantiate a filterset and save it as an attribute
        # on the view instance for later.
        # Return the filtered queryset
        order=self.get_order_by()
        self.filter_count = self.filterset.qs.distinct().count()
        if order:           
            order = order.replace(".", "__")
            return self.filterset.qs.distinct().order_by(order)
    
        return self.filterset.qs.distinct()

    #--------------------------------------------------------------------------
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        self.context_list = context['object_list']
                
        filter_record_dict = self.get_filter_record()
        filter_record = "Selected: "+ str(filter_record_dict).replace("{", "").replace("}", "") if str(filter_record_dict).replace("{", "").replace("}", "") else None

        # Pass the filterset to the template - it provides the form.
        # self.filterset.update_choice_filters(filter_record_dict)
        
        context['filter'] = self.filterset
        context['paginate_by'] = self.get_paginate_by(self, **kwargs)
        context['fields'] = self.model.get_fields(fields = self.model_fields)
        context['filterset'] = filter_record
        context['Count'] = self.model.objects.count()
        context['querycount'] = self.filter_count

        return context
    
    #--------------------------------------------------------------------------
    def get_paginate_by(self, queryset):
        qs=super().get_queryset()
        paginate_by= self.request.GET.get("paginate_by", self.paginate_by)

        return paginate_by
 
    #--------------------------------------------------------------------------
    def get_order_by(self):   

        order_by=self.request.GET.get("order_by", self.order_by) or None
        acs_decs=""
        if order_by:
            order_field=""
            if order_by[0]=="-":
                acs_decs=order_by[0]
                order_field=order_by[1:]
            else:
                order_field=order_by              
            index=self.find_item_index(list(self.model_fields.values()), order_field)
            order_by=acs_decs+ list(self.model_fields.keys())[index]
            return order_by
        
        return order_by 
    

# -----------------------------------------------------------------
# -- View for simple update files and images to database--
# -----------------------------------------------------------------
class File_CreateView(LoginRequiredMixin,FormView):
    form_class = None
    model = None
    file_field = None
    related_name = None
    transaction_use = 'default'
    transaction_use_manytomany = 'default'

    def dispatch(self, request, *args, **kwargs):
        self.object_ = get_object_or_404(self.model, organism_id=kwargs['pk'])
        return super().dispatch(request, *args, **kwargs)

    def post(self, request, *args, **kwargs):
        # Handle AJAX file upload
        if request.headers.get('x-requested-with') == 'XMLHttpRequest':
            file_data = request.FILES.get(self.file_field)
            if file_data:           
                file_name =file_data.name
                file_type = file_data.content_type
                
                response = {
                    'name': file_name,
                    'type': file_type,
                }
                return JsonResponse(response)
        # Handle form submission
        else:
            return super().post(request, *args, **kwargs)

    def form_valid(self, form):
        instance = form.save(commit=False)
        if getattr(instance, self.file_field):
            with transaction.atomic(using=self.transaction_use):
                kwargs={'user': self.request.user}
                instance.save(**kwargs)
            with transaction.atomic(using=self.transaction_use_manytomany):
                getattr(self.object_, self.related_name).add(instance)
                self.object_.save(**kwargs)            
        else:
            messages.warning(self.request, 'No file provided.')
        return redirect(self.request.META['HTTP_REFERER'])

    def form_invalid(self, form):
        messages.warning(self.request, f'Update failed due to {form.errors} error')
        return redirect(self.request.META['HTTP_REFERER'])

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["object"] = self.object_
        return context
    
# -----------------------------------------------------------------
# -- Export View class--
# -----------------------------------------------------------------

class Base_DataExportView(LoginRequiredMixin, View):
    '''
    Export Data in EXCEL or CSV Format
    '''
    selected_pks_string = None
    organism = None    # organism PK value - data from a table in organism detail view

    def post(self, request):
        items=None
        app_name = request.POST.get('app_name')
        model_name = request.POST.get('model_name')
        current_date = datetime.now().strftime('%Y-%m-%d')
        filename = f'{model_name}_{current_date}'
        
        if not app_name or not model_name:
            messages.error(request, "No application or model name were provided for export.")
            return redirect(request.path)
        try:
            Model = apps.get_model(app_name, model_name)
        except LookupError:
            # Handle the case where the model does not exist.
            return HttpResponse("Model not found.")
        
        self.selected_pks_string = request.POST.get('selected_pks')
        self.organism = request.POST.get('organism_pk')

        if self.selected_pks_string == 'SelectAll':
            items = Model.objects.filter(pk__in=request.session.get(f"{Model}_cached_queryset") or Model.objects.all())
            # table = Model.get_pivottable(querydata=items, columns_str=columns_str, index_str=index_str, aggfunc=aggfunc_name, values=values)
        elif self.selected_pks_string:
            selected_pks = json.loads(self.selected_pks_string)
            items = Model.objects.filter(pk__in=selected_pks)
        else:
            if self.organism:
                displaycols = ['Drug Class', 'Drug Name', 'MIC', 'BP Profile', 'BatchID', 'Source', 'BP Source']
                df = get_Antibiogram_byOrgID(self.organism)
                pivdf = piv_Antibiogram_byOrgID(df)
                df.reset_index(inplace=True)
                df = df[displaycols]
                filename = f'Antibiogram_{current_date}'
            else:
                messages.error(request, "No items were selected for export.")
                return redirect(request.META['HTTP_REFERER'])
        if items:
            df = pd.DataFrame.from_records(items.values())

        if 'acreated_at' in df.columns:
            df['acreated_at'] = pd.to_datetime(df['acreated_at']).dt.tz_localize(None).apply(lambda a: a.date())
        if 'aupdated_at' in df.columns:
            df['aupdated_at'] = pd.to_datetime(df['aupdated_at']).dt.tz_localize(None).apply(lambda a: a.date())
        if 'adeleted_at' in df.columns:
            df['adeleted_at'] = pd.to_datetime(df['adeleted_at']).dt.tz_localize(None).apply(lambda a: a.date())

        response = HttpResponse("No valid export option was selected.")
        if "csvdownload" in request.POST:
            response = HttpResponse(content_type='text/csv')
            response['Content-Disposition'] = f'attachment; filename="{filename}.csv"'
            df.to_csv(path_or_buf=response, index=False)

        elif "xlsxdownload" in request.POST:
            response = HttpResponse(content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
            response['Content-Disposition'] = f'attachment; filename="{filename}.xlsx"'
            df.to_excel(excel_writer=response, index=False)

        return response