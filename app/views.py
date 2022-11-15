from django.shortcuts import render, redirect, get_object_or_404
from .models import ApplicationUser, Dictionaries
from aa_chem.models import Drugbank, Taxonomy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.contrib import messages
from django.utils.decorators import method_decorator
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.urls import reverse_lazy, reverse
from .forms import ApplicationUser_form, Dictionary_form
from django.contrib.auth import logout
from django.contrib.auth.models import Permission
# Create your views here.
# def check_admin(user):
#    return user.is_superuser

# @user_passes_test(check_admin)
# def my_view(request): 


# =================================Login into Home Page==================================
def index(req):
    if req.user.is_authenticated:

        user=ApplicationUser.objects.get(username=req.user.username)
        # print(user.username)
        
        if user.is_appuser==False:
            logout(req)
            user.delete()
            messages.warning(req, 'User not authorized, please contact Admin!')
            return redirect("/")
    return render(req, 'aa_chem/home.html')



# =========================Application Users View========================================
class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):

    def test_func(self):
        return self.request.user.is_superuser

# @permission_required('app.change_groupfilter')
@login_required(login_url='/login/')
def userprofile(req, id):
    current_user=get_object_or_404(User, pk=id)
    permissions = Permission.objects.filter(user=current_user)
    return render(req, 'app/userprofile.html', {'currentUser': current_user, 'perm':permissions})


class AppUserListView(LoginRequiredMixin, ListView):
    model=ApplicationUser
    fields='__all__'
    template_name = 'app/appUsers.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
       
        print(context["objects"])
       
        return context

class AppUserCreateView(SuperUserRequiredMixin, CreateView):
    model=ApplicationUser
    fields=['user_id', 'permissions', 'is_appuser']
    template_name = 'app/appUsersCreate.html'
    success_url = reverse_lazy('userslist')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context


class AppUserUpdateView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    fields=['user_id', 'permissions', 'is_appuser']
    template_name = 'app/appUsersUpdate.html'
    success_url = reverse_lazy('userslist')

    def form_valid(self, form):

        if form.instance.permissions=='staff':
            form.instance.is_staff=True
            form.instance.is_superuser=False
            
        if form.instance.permissions=='admin':
            form.instance.is_staff=True
            form.instance.is_superuser=True
        
        return super().form_valid(form)

class AppUserDeleteView(SuperUserRequiredMixin, DeleteView):
    model=ApplicationUser
    template_name='app/appUsersDel.html'
    success_url = reverse_lazy('userslist')

# =====================================Dictionaries View==========================

class DictionariesView(LoginRequiredMixin, ListView):
    model=Dictionaries
    fields="__all__"
    template_name='app/dictList.html'
    
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
        # else:
            form_error=True    
            messages.error(req, form.errors)
            return redirect("dict_view")
    
    return render(req, 'app/dictCreate.html', {'form': form, 'form_error':form_error})

@user_passes_test(lambda u: u.is_superuser, redirect_field_name=None) #login_url='/redirect/to/somewhere'
def deleteDictionary(req, pk):
    kwargs={}
    kwargs['user']=req.user
    context={}
    object_=get_object_or_404(Dictionaries, Dict_Value=pk)
    context['object']=object_
    if req.method=="POST":      
        object_.delete(**kwargs)
        print("deleted")
        return redirect(req.META['HTTP_REFERER'])
    
    return render(req, 'app/dictionary_del.html', context)