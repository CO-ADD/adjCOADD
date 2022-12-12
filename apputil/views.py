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
from dorganism.utils import import_excel
from dorganism.models import Organism, Taxonomy

# ==========utilized in Decoration has_permissions, an Alert on Permissions ==========
def permission_not_granted(req):
    return HttpResponse("Permission Not Granted")



## =================================APP Home========================================
def index(req):
    object_1=Organism.objects.count()
    object_2=Taxonomy.objects.count()
    return render(req, 'dorganism/home.html', {'objects_org': object_1, 'objects_taxo':object_2})
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
                return redirect("index")
        else:
            form = Login_form()
        return render(req, 'registration/login.html', {'form': form})    

def logout_user(req):
    logout(req)    
    return redirect("/")

# =================================APP Log in/out ==================================##

## =========================Application Users View====================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):

    def test_func(self):
        return self.request.user.is_superuser

@login_required(login_url='/login/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    return render(req, 'apputil/userprofile.html', {'currentUser': current_user})

class AppUserListView(LoginRequiredMixin, ListView):
    model=ApplicationUser
    fields='__all__'
    template_name = 'apputil/appUsers.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context

class AppUserCreateView(SuperUserRequiredMixin, CreateView):
    model=ApplicationUser
    fields=['name', 'username', 'permission']
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context

class AppUserUpdateView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    fields=['name', 'permission', ]
    template_name = 'apputil/appUsersUpdate.html'
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

class DictionaryView(LoginRequiredMixin, ListView):
    model=Dictionary
    fields="__all__"
    template_name='apputil/dictList.html'
    
    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.filter(astatus__gte=0)
        return context
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



# ==========================File Process=====================================

from django.core.files.storage import FileSystemStorage
from django.views import View
from django import forms
import json
from django.core import serializers
import os

class FileUploadForm(forms.Form):
    excel_file = forms.FileField(required=False)

class Importhandler(View):
    template_name='apputil/importdata.html'
    file_url=''
    data_list=[]

    def get(self, request):
        form = FileUploadForm()
        return render(request, 'apputil/importdata.html', { 'form': form, })

    def post(self, request):
        form = FileUploadForm(request.POST, request.FILES)
        context = {}

        try:
            if form.is_valid():
                myfile=request.FILES['myfile']
                fs=FileSystemStorage()
                filename=fs.save(myfile.name, myfile)
                self.file_url=fs.url(filename)
                context['message']=self.file_url
         
            
                return render(request,'apputil/importdata.html', context)
            else:
                messages.warning(request, f'There is {form.errors} error, upload again')
        except Exception as err:
            messages.warning(request, f'There is {err} error, upload again. myfile error-- filepath cannot be null, choose a correct file')

        context['form'] = form
       

        return render(request, 'apputil/importdata.html', context)
    
    #delete task
    @csrf_exempt
    @staticmethod
    def delete_task(request):
        print(Importhandler.file_url)
        for filename in os.listdir("uploads"):
                file_path=os.path.join("uploads", filename)
                try:
                    os.unlink(file_path)
                    print("removed!")
                except Exception as err:
                    print(err)
        return JsonResponse({})

    @csrf_exempt
    @staticmethod
    def run_task(request):
        
        task_type = request.POST.get("type")
        Importhandler.file_url=request.POST.get("filepath")
        
        if task_type =='Validation':
            # read file from self.file_url
            # validate fields
            #return status=error, warning, or pass
            #return result=ErrorList
            
            Importhandler.data_list=import_excel(Importhandler.file_url)
            if isinstance(Importhandler.data_list, Exception):
                result=str(Importhandler.data_list)
                status='Form Errors'
            else:
                result='pass'
                status='Form is Valid'
            return JsonResponse({"task_num": "somehash", 'task_status':status, 'task_result': result})
        elif task_type=='Cancel':
            Importhandler.delete_task(request)
            return JsonResponse({})

    def proceed_save(request):
        for obj in Importhandler.data_list:
            obj.save()
        return JsonResponse({"status":"SUCCESS"}, status=200)
        # pass


    
    # create entries
    @csrf_exempt
    @staticmethod
    def save_task(request):
        # call save object funtion 
        pass
    
