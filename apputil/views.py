from django.contrib.admin.views.decorators import staff_member_required
from django.contrib import messages
from django.contrib.auth import logout, login
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.contrib.auth.forms import AuthenticationForm
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.shortcuts import get_object_or_404, HttpResponse, render, redirect
from django.http import JsonResponse
from django.urls import reverse_lazy, reverse
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.views.generic import ListView, TemplateView
from django.views.generic.edit import UpdateView, CreateView, DeleteView

from .forms import ApplicationUser_form, Dictionary_form, Login_form
from .models import ApplicationUser, Dictionary
from .utils_dataimport import import_excel
from dorganism.models import Organism, Taxonomy
from apputil.utils import FilteredListView
from dorganism.utils import Dictionaryfilter
from adjcoadd.constants import *
from ddrug.models import VITEK_Card, VITEK_ID, VITEK_AST
from .utils import SuperUserRequiredMixin

# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")




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
    return render(req, 'apputil/userprofile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, ListView):
    login_url = '/'
    model=ApplicationUser
    fields='__all__'
    template_name = 'apputil/appUsers.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context
        
from django.http import QueryDict
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
    fields=['name', 'username', 'permission',]
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')


class AppUserDeleteView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    template_name='apputil/appUsersDel.html'
    success_url = reverse_lazy('userslist')
    fields=['is_appuser']

    def form_valid(self, form):

        form.instance.is_appuser==False
        print(form.instance.is_appuser)
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
@user_passes_test(lambda u: u.has_permission('Admin'), redirect_field_name=None)
def createDictionary(req):
    kwargs={}
    kwargs['user']=req.user
    form=Dictionary_form()
    form_error=False
    if req.method=='POST':
        form=Dictionary_form(req.POST)
        try:
            if form.is_valid:
                print("form is valid")
                instance=form.save(commit=False)
                instance.save(**kwargs)
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
    if req.headers.get('x-requested-with') == 'XMLHttpRequest' and req.method == "POST":
        dict_value=req.POST.get("dict_value")
        dict_class=req.POST.get("dict_class")
        dict_desc=req.POST.get("dict_desc")
        object_=get_object_or_404(Dictionary, dict_value=dict_value)
        try:
            if object_:
                object_.dict_class=dict_class
                object_.dict_desc=dict_desc
                print(object_.dict_desc)
                object_.save(**kwargs)
                return JsonResponse({"result": "Saved"})
        except Exception as err:
            return JsonResponse({"result": err})
    
    return JsonResponse({})

