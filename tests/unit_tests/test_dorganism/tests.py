import threading
from django.test import RequestFactory, TestCase, Client, TransactionTestCase
from django.urls import reverse
# from django_multitenant.utils import set_current_schema
from dorganism.models import Organism, Taxonomy
from dorganism.views import *



# ------------------
#Test Concurrency when Create Entry
# -----------------
from multiprocessing import Process
from django.contrib.auth.models import User
from apputil.models import ApplicationUser, Dictionary


# class OrganismCreationTestCase(TestCase):

#     def setUp(self):
#         self.client = Client()
#         self.user = ApplicationUser.objects.create_user(username='test', name='test', permission='Admin')
#         self.client.login(username='test', password='Django')

    # def test_concurrent_creation(self):
    #     url = reverse('org_create')
    #     data = {'search_organism': 'Test Organism'}
        
    #     def create_organism():
    #         self.client.post(url, data)

    #     # Start two concurrent requests to create an organism
    #     p1 = Process(target=create_organism)
    #     p2 = Process(target=create_organism)
    #     p1.start()
    #     p2.start()
    #     p1.join()
    #     p2.join()

    #     # Check that only one organism was created
    #     count = Organism.objects.filter(name='Test Organism').count()
    #     self.assertEqual(count, 1)


# def test_concurrently(times):
#     """
#     Add this decorator to small pieces of code that you want to test
#     concurrently to make sure they don't raise exceptions when run at the
#     same time.  E.g., some Django views that do a SELECT and then a subsequent
#     INSERT might fail when the INSERT assumes that the data has not changed
#     since the SELECT.
#     """
#     def test_concurrently_decorator(test_func):
#         def wrapper(*args, **kwargs):
#             exceptions = []
#             def call_test_func():
#                 try:
#                     test_func(*args, **kwargs)
#                 except Exception as e:
#                     exceptions.append(e)
#                     raise
#             threads = []
#             for i in range(times):
#                 threads.append(threading.Thread(target=call_test_func))
#             for t in threads:
#                 t.start()
#             for t in threads:
#                 t.join()
#             if exceptions:
#                 raise Exception('test_concurrently intercepted %s exceptions: %s' % (len(exceptions), exceptions))
#         return wrapper
#     return test_concurrently_decorator


#
# Test Create View
TransactionTestCase.databases = {"default", "survey"}
class CreateTaxonomyTestCase(TransactionTestCase):
    # multi_db = True
    # databases={'dorganism', 'default'}

    def setUp(self):
        self.factory=RequestFactory()
        self.user= ApplicationUser.objects.create(name='orgdb', username='orgdb', permission='Admin')
        self.organism=Taxonomy.objects.create(organism_name='for_test', urlname='fortest', tax_id=1, parent_tax_id=2000, tax_rank='test', code='tetcode', other_names='test' )
        # self.dictionary=Dictionary.objects.create(dict_value='test_unique', dict_class='org_class', dict_app='test', dict_desc='testcase')
        
        self.create_url=reverse('taxo_create')   
    
    def test_create(self):
        request=self.factory.get(self.create_url, follow=True)
        request.user = self.user
        request.organism=self.organism
        response = TaxonomyCreateView.as_view()(request)

        self.assertEqual(response.status_code, 200)
    # def test_usercreate_GET(self):
    #     request = self.factory.get(self.create_url, follow=True)
    #     request.user = self.user
    #     response=AppUserCreateView.as_view()(request)
    #     self.assertEquals(response.status_code, 200)