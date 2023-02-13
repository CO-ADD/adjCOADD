import json
import time
from django.test import RequestFactory,TestCase, Client
from django.urls import reverse
from apputil.models import *
from django.contrib.auth.models import AnonymousUser
from apputil.views import *

TestCase.databases={'default'}

class TestAppUserView(TestCase):
       
    def setUp(self):
        self.factory = RequestFactory()
        self.user= ApplicationUser.objects.create_user(name='testcase1', username='testcase1', permission='Admin')
    
        self.client=Client()
        self.list_url=reverse('userslist')
        self.create_url=reverse('createAppUser')   
        self.update_url=reverse('updateAppUser', args=[self.user.pk])
        self.delete_url=reverse('deleteAppUser', args=[self.user.pk])
    
    def test_userlist_GET(self):
        response=self.client.get(self.list_url, follow=True)
       
        print(response)
        self.assertEquals(response.status_code, 200)

    def test_usercreate_GET(self):
        request = self.factory.get(self.create_url, follow=True)
        request.user = self.user
        response=AppUserCreateView.as_view()(request)
        self.assertEquals(response.status_code, 200)
    
    def test_userupdate_GET(self):
        request = self.factory.get(self.update_url, follow=True)
        request.user = self.user
        response=updateApplicationuser(request, request.user.pk)
        self.assertEquals(response.status_code, 200)

    def test_userdelete_GET(self):
        request = self.factory.get(self.delete_url, follow=True)
        request.user = self.user
        request.user.pk=self.user.pk
        response=AppUserDeleteView.as_view()(request, pk=request.user.pk)
        self.assertEquals(response.status_code, 200)