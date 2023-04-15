from django.test import SimpleTestCase
from django.urls import reverse, resolve

from apputil.views import *

class TestUrls(SimpleTestCase):

    def test_login_url_is_resolved(self):
        url=reverse('login')
        self.assertEquals(resolve(url).func, login_user)
    
    def test_AppUserList_url_is_resolved(self):
        url=reverse('userslist')
        self.assertEquals(resolve(url).func.view_class, AppUserListView)

    def test_AppUserUpdate_url_is_resolved(self):
        url=reverse('updateAppUser', args=['somePk'] )
        self.assertEquals(resolve(url).func, updateApplicationuser)