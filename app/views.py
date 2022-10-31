from django.shortcuts import render, redirect, get_object_or_404
from .models import ApplicationUser, Dictionaries
from aa_chem.models import Drugbank, Taxonomy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.utils.decorators import method_decorator
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.urls import reverse_lazy, reverse
from .forms import CreateApplicationUser_form, Dictionary_form
from django.contrib.auth import logout
from django.contrib.auth.models import Permission
# Create your views here.
# def check_admin(user):
#    return user.is_superuser

# @user_passes_test(check_admin)
# def my_view(request): 

def index(req):
    
    objects= Taxonomy.objects.all()
    if req.user.is_authenticated:

        user=ApplicationUser.objects.get(username=req.user.username)
        print(user.username)
        
        if user.is_appuser==False:
            logout(req)
            user.delete()
            return redirect("/")
    return render(req, 'aa_chem/home.html', {'objects': objects})

class SuperUserRequiredMixin(LoginRequiredMixin, UserPassesTestMixin):

    def test_func(self):
        return self.request.user.is_superuser

# @permission_required('app.change_groupfilter')
@login_required
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
    # fields='__all__'
    form_class=CreateApplicationUser_form
    template_name = 'app/appUsersCreate.html'
    success_url = reverse_lazy('userslist')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context


class AppUserUpdateView(SuperUserRequiredMixin, UpdateView):
    model=ApplicationUser
    fields='__all__'
    template_name = 'app/appUsersUpdate.html'
    success_url = reverse_lazy('userslist')

class AppUserDeleteView(SuperUserRequiredMixin, DeleteView):
    model=ApplicationUser
    template_name='app/appUsersDel.html'
    success_url = reverse_lazy('userslist')



class DictionariesView(ListView):
    model=Dictionaries
    fields="__all__"
    template_name='app/dictList.html'
    
    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
       
        print(context["objects"])
       
        return context

def DictCreate(req):
    form=Dictionary_form()
    if req.method=='POST':
        form=Dictionary_form(req.POST)
        if form.is_valid:
            instance=form.save(commit=False)
            instance.save(req.user)
            return redirect("dict_view")
    return render(req, 'app/dictCreate.html', {'form': form})