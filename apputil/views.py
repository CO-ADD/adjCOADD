from django.shortcuts import render, redirect, get_object_or_404
from .models import ApplicationUser, Dictionary
from dorganism.models import  Taxonomy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.contrib import messages
from django.utils.decorators import method_decorator
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.urls import reverse_lazy, reverse
from .forms import ApplicationUser_form, Dictionary_form, Login_form
from django.contrib.auth import logout, login
from django.contrib.auth.models import Permission
from django.contrib.auth.forms import AuthenticationForm
# Create your views here.
# def check_admin(user):
#    return user.is_superuser

# @user_passes_test(check_admin)
# def my_view(request): 


# =================================Login into Home Page==================================
def login_user(req):
    if req.method=='POST':
        form=Login_form(data=req.POST)
        username_ldap=req.POST.get('username')
        # print(user)
        if form.is_valid():
            print("form is valid")
            user=form.get_user()
            login(req, user, backend="django_auth_ldap.backend.LDAPBackend",)
            return redirect("/")
        
        else:
            messages.warning(req, ' no permission for this application, please contact Admin!')
            return redirect(req.META['HTTP_REFERER'])
    else:
        form = Login_form()
    return render(req, 'registration/login.html', {'form': form})    

def logout_user(req):
    logout(req)    
    return redirect("/accounts/login/")



def index(req):
    return render(req, 'dorganism/home.html')



# =========================Application Users View========================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):

    def test_func(self):
        return self.request.user.is_superuser


@login_required(login_url='/login/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    permissions = Permission.objects.filter(user=current_user)
    return render(req, 'apputil/userprofile.html', {'currentUser': current_user, 'perm':permissions})


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
    fields=['user_id', 'permissions']
    template_name = 'apputil/appUsersCreate.html'
    success_url = reverse_lazy('userslist')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context


class AppUserUpdateView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    fields=['user_id', 'permissions', ]
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


# =====================================Dictionary View==========================

class DictionaryView(LoginRequiredMixin, ListView):
    model=Dictionary
    fields="__all__"
    template_name='apputil/dictList.html'
    
    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.filter(astatus__gte=0)
        return context
# 
@user_passes_test(lambda u: u.is_superuser, redirect_field_name=None)
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

# @user_passes_test(lambda u: u.is_superuser, redirect_field_name=None) #login_url='/redirect/to/somewhere'
# def deleteDictionary(req, pk):
#     kwargs={}
#     kwargs['user']=req.user
#     context={}
#     object_=get_object_or_404(Dictionary, Dict_Value=pk)
#     context['object']=object_
#     if req.method=="POST":      
#         object_.delete(**kwargs)
#         print("deleted")
#         return redirect(req.META['HTTP_REFERER'])
    
#     return render(req, 'apputil/dictionary_del.html', context)