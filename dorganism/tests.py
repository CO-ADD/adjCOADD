import threading
from django.test import RequestFactory, TestCase
from .models import Organism, Taxonomy
from .views import updateOrganism, createOrgnisms, TaxoCreate

TestCase.databases={'dorganism'}
# ====================================test Concurrently=======================================
def test_concurrently(times):
    """
    Add this decorator to small pieces of code that you want to test
    concurrently to make sure they don't raise exceptions when run at the
    same time.  E.g., some Django views that do a SELECT and then a subsequent
    INSERT might fail when the INSERT assumes that the data has not changed
    since the SELECT.
    """
    def test_concurrently_decorator(test_func):
        def wrapper(*args, **kwargs):
            exceptions = []
            def call_test_func():
                try:
                    test_func(*args, **kwargs)
                except Exception as e:
                    exceptions.append(e)
                    raise
            threads = []
            for i in range(times):
                threads.append(threading.Thread(target=call_test_func))
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            if exceptions:
                raise Exception('test_concurrently intercepted %s exceptions: %s' % (len(exceptions), exceptions))
        return wrapper
    return test_concurrently_decorator


# class CreateOrganismTestCase(TestCase):
#     def setUp(self):
#         self.factory=RequestFactory()
#         self.organism=Organism.objects.create_orgnisms(Organism_Name='for test', Risk_Group='RG03' )

#     @test_concurrently(15)
#     def test_update(self):
#         request=self.factory.get('/organism/createOrg/')
#         request.organism=self.organism
#         response = createOrgnisms(request)

#         self.assertEqual(response.status_code, 200)


# Test Create View

def test_Organism_create(self):
    self.client.post('createOrg/', {'author':"manualvarado22", 'title': "Super Important Test", 'content':"This is really important.", 'published_date':timezone.now()})
    self.assertEqual(Post.objects.last().title, "Super Important Test")

def test_display_post(self):
    post = Post.objects.create(...whatever...)
    response = self.client.get(reverse('blog:post_detail', pk=post.pk))
    self.assertContains(response, "really important")

class CreateTaxonomyTestCase(TestCase):
    def setUp(self):
        self.factory=RequestFactory()
        self.organism=Taxonomy.objects.create_taxonomy(organism_name='for test', risk_group='RG03', lineaage=['Lineage', 'test'], tax_id=1, parent_tax_id=2000, tax_rank='test', code='tetcode', other_names='test' )

    @test_concurrently(15)
    def test_update(self):
        request=self.factory.get('/organism/createOrg/')
        request.organism=self.organism
        response = TaxoCreate(request)

        self.assertEqual(response.status_code, 200)

