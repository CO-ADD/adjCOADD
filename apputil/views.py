from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse, QueryDict
from django.http import JsonResponse, QueryDict
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView, CreateView, DeleteView

from adjcoadd.constants import *
from dorganism.models import Organism, Taxonomy

from .forms import ApplicationUser_form, Dictionary_form, Login_form
from .models import ApplicationUser, Dictionary
from .utils import SuperUserRequiredMixin, permission_not_granted, FilteredListView, AppUserfilter, Dictionaryfilter





## =================================APP Home========================================

# import setup
@login_required(login_url='/')
def index(req):
    # print(setup.version)
    object_1=Organism.objects.count()
    object_2=Taxonomy.objects.count()
    
    return render(req, 'dorganism/home.html', {'objects_org': object_1, 'objects_taxo':object_2,})
## =================================APP Home======================================##

## =================================APP Log in/out =================================
def login_user(req):
    if req.user.is_authenticated:
        return redirect("index")
    else:
        if req.method=='POST':
            form=Login_form(data=req.POST)
            username_ldap=req.POST.get('username')
        # print(user)
            if form.is_valid():
                print("form is valid")
                user=form.get_user()
                login(req, user, backend="django_auth_ldap.backend.LDAPBackend",)
                return redirect("index")
        
            else:
                messages.warning(req, ' no permission for this application, please contact Admin!')
                return redirect("/")
        else:
            form = Login_form()
        return render(req, 'registration/login.html', {'form': form})    

def logout_user(req):
    logout(req)    
    return redirect("/")

# =================================APP Log in/out ==================================##

## =========================Application Users View====================================

@login_required(login_url='/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    return render(req, 'apputil/appUserProfile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=ApplicationUser
    template_name = 'apputil/appUsers.html'  
    filterset_class = AppUserfilter
    model_fields=APPUSER_FIELDs


@user_passes_test(lambda u: u.has_permission('Admin'), login_url='permission_not_granted') 
def updateApplicationuser(req, pk):
    object_=get_object_or_404(ApplicationUser, pk=pk)
    form=ApplicationUser_form(instance=object_)
    context={
        "form":form,
        "object":object_,
        }
    if req.method=='PUT':
        qd=QueryDict(req.body).dict()       
        form=ApplicationUser_form(data=qd, instance=object_)
        if form.is_valid():
            instance=form.save()
            context={"object":object_}
            return render(req, "apputil/appuser_tr.html", context)
    return render(req, "apputil/appUsersUpdate.html", context)

class AppUserCreateView(SuperUserRequiredMixin, CreateView):
    model=ApplicationUser
    form_class = ApplicationUser_form
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')


class AppUserDeleteView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    template_name='apputil/appUsersDel.html'
    success_url = reverse_lazy('userslist')
    fields=['is_appuser']

    def form_valid(self, form):
        form.instance.is_appuser==False
        return super().form_valid(form)
# =========================Application Users View==================================##

## ========================Dictionary View===========================================

class DictionaryView(LoginRequiredMixin, FilteredListView):
    login_url = '/'
    model=Dictionary
    template_name='apputil/dictList.html'
    filterset_class = Dictionaryfilter
    model_fields=DICTIONARY_FIELDs

    
# 
@user_passes_test(lambda u: u.has_permission('Admin'), login_url='/')
def createDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    form=Dictionary_form()
    form_error=False
    if req.method=='POST':
        form=Dictionary_form(req.POST)
        try:
            if form.is_valid:
                instance=form.save(commit=False)
                instance.save(**kwargs)
                print("save")
                return redirect("dict_view")
        except Exception as err:
       
            form_error=True    
            messages.error(req, form.errors)
            return redirect("dict_view")
    
    return render(req, 'apputil/dictCreate.html', {'form': form, 'form_error':form_error})
## ============================Dictionary View======================================##
@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def updateDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    if req.user.has_permission('Admin'):
        print(req.user.has_permission('Admin'))

        if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
            dict_value=req.POST.get("dict_value") 
            type=req.POST.get("type") or None
            print(type)
            value=req.POST.get("value").strip() or None
            print(value)
            object_=get_object_or_404(Dictionary, dict_value=dict_value)
            try:
                if object_:
                    if type=='dict_class':
                        object_.dict_class=value
                    if type=='dict_desc':
                        object_.dict_desc=value
                    print(object_.dict_desc)
                    object_.save(**kwargs)
                    return JsonResponse({"result": "Saved"})
            except Exception as err:
                return JsonResponse({"result": err})
    else:
        return JsonResponse({"result": "Permission_denied"})
    return JsonResponse({})


@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def deleteDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
        dict_value=req.POST.get("dict_value") 
        object_=get_object_or_404(Dictionary, dict_value=dict_value)
        try:
            if object_:
                object_.delete(**kwargs)
                return JsonResponse({"result": "Deleted!"})
        except Exception as err:
            return JsonResponse({"result": err})
    
    return JsonResponse({})


############################################### Export CSV View ###########################################
import csv
import datetime 
from django.apps import apps

@login_required
@user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
def exportCSV(request):
    if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
        data_arr = request.POST.getlist('data_arr[]')
        data_fields = request.POST.getlist('fields[]')
        model_name=request.POST.get('model_name')
        try:
            model=apps.get_model('dorganism', model_name)
        except:
            model=apps.get_model('ddrug', model_name)
        query=model.objects.filter(pk__in=data_arr)
        response = HttpResponse(content_type='text/csv')
        file_name = "fltred_loaction_data" + str(datetime.date.today()) + ".csv"
        writer = csv.writer(response)
        writer.writerow(data_fields)
        for i in query.values_list(*data_fields):
            writer.writerow(i)
        response['Content-Disposition'] = 'attachment; filename = "' + file_name + '"'
        return response


############################################### Export CSV View ###########################################
import csv
import datetime 
from django.apps import apps

@login_required
@user_passes_test(lambda u:u.has_permission('Admin'), login_url='permission_not_granted') 
def exportCSV(request):
    if request.headers.get('x-requested-with') == 'XMLHttpRequest' and request.method == "POST":
        data_arr = request.POST.getlist('data_arr[]')
        data_fields = request.POST.getlist('fields[]')
        model_name=request.POST.get('model_name')
        try:
            model=apps.get_model('dorganism', model_name)
        except:
            model=apps.get_model('ddrug', model_name)
        query=model.objects.filter(pk__in=data_arr)
        response = HttpResponse(content_type='text/csv')
        file_name = "fltred_loaction_data" + str(datetime.date.today()) + ".csv"
        writer = csv.writer(response)
        writer.writerow(data_fields)
        for i in query.values_list(*data_fields):
            writer.writerow(i)
        response['Content-Disposition'] = 'attachment; filename = "' + file_name + '"'
        return response

# ==============Test Site: Combined Queries of Models ==========================#
from django.core.paginator import Paginator
# from itertools import chain

def testsite(request):
    choice_class=Dictionary.objects.order_by().values('dict_class').distinct()
    a=[tuple(d.values()) for d in choice_class]
    b=[(x[0], x[0]) for x in a]
    print(b)
    context={}
    # organism = Organism.objects.values_list('organism_name')
    # taxonomy = Taxonomy.objects.values_list('organism_name')
    objects=[]#organism.union(taxonomy).order_by('-pk')
    if request.method == 'POST':
        searched = request.POST['searched']
        organism = Organism.objects.filter(organism_name__organism_name__icontains=searched,
                                   ).values_list("organism_name")
        taxonomy = Taxonomy.objects.filter(organism_name__icontains=searched,
                                      ).values_list('organism_name')
       
        # objects= [item for item in chain(organism, taxonomy)]
        objects=organism.union(taxonomy).order_by('-pk')
        print(objects)

    paginator = Paginator(objects, 2) # Show 25 contacts per page.

    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)

    print(page_obj)

    context = {
        'page_obj': page_obj,
    }

        
    
    return render(request, "utils/test.html", context)
 