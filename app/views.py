from django.shortcuts import render, redirect, get_object_or_404
from .models import User, Groupfilter,ApplicationUser
from aa_chem.models import Drugbank, Taxonomy
from django.views.generic.edit import UpdateView, CreateView, DeleteView
from django.views.generic import ListView
from django.utils.decorators import method_decorator
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.contrib.admin.views.decorators import staff_member_required
from django.contrib.auth.decorators import user_passes_test, login_required, permission_required
from django.urls import reverse_lazy
from .forms import GroupCreate 
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

        user=User.objects.get(username=req.user.username)
        
        if user.role=="":
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


class UserListView(LoginRequiredMixin, ListView):
    model=ApplicationUser
    fields='__all__'
    template_name = 'app/appUsers.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
       
        print(context["objects"])
       
        return context

class GroupListView(LoginRequiredMixin, ListView):
    model=Groupfilter
    fields='__all__'
    template_name = 'app/groups.html'

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        context["objects2"]=User.objects.all()
        print(context["objects"])
        print(context["objects2"])
        return context

class GroupCreateView(SuperUserRequiredMixin, CreateView):
    model=Groupfilter
    # fields='__all__'
    form_class=GroupCreate
    template_name = 'app/groupsCreate.html'
    success_url = reverse_lazy('usermanage')

    def get_context_data(self, **kwargs):
        context=super().get_context_data(**kwargs)
        context["objects"]=self.model.objects.all()
        return context


class GroupUpdateView(SuperUserRequiredMixin, UpdateView):
    model=Groupfilter
    fields='__all__'
    template_name = 'app/groupsUpdate.html'
    success_url = reverse_lazy('usermanage')

class GroupDeleteView(SuperUserRequiredMixin, DeleteView):
    model=Groupfilter
    template_name='app/groupsDelete.html'
    success_url = reverse_lazy('usermanage')