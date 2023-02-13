import json
import time
from django.test import TestCase, Client
from django.urls import reverse
from apputil.models import *


TestCase.databases={'default'}

class TestAppUserView(TestCase):
    # @classmethod
    # def setUpTestData(cls):
    #     # Set up data for the whole TestCase
    #     cls.appusers = ApplicationUser.objects.create(name="Test"+str(time.time()), username="Test"+str(time.time()), permission="Admin")
    
    def setUp(self):
        self.client=Client()
        self.list_url=reverse('userslist')
        self.update_url=reverse('updateAppUser', args=['test'])
        # morelinksURL
        ApplicationUser.objects.create(name='testcase1', username='testcase1', permission='Admin')
    
    def test_userlist_GET(self):
        response=self.client.get(self.list_url, follow=True)
       
        print(response)
        self.assertEquals(response.status_code, 200)
        # self.assertTemplateUsed(response, 'apputil/appUsers.html')
    
    # def test_userlist_GET(self):
    #     response=self.client.get(self.update_url)
       
    #     print(response)