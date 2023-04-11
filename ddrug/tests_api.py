from django.test import TestCase
from apputil.models import ApplicationUser as User
import requests
from django.urls import reverse

class TestApiAccess(TestCase):
    # databases={'ddrug', 'default'}
    def setUp(self):
        self.username = 'uqwzhon4'
        self.password = 'Tddx^_^10'

    def tearDown(self):
        pass

    def test_obtain_access_token(self):
        url = reverse('token_obtain_pair')
        credentials = {
            'username': self.username,
            'password': self.password,
        }

        response = self.client.post(url, data=credentials)
        if 'access' not in response.json():
            print("Failed to obtain access token:", response.json())
        self.assertEqual(response.status_code, 200)
        self.assertIn('access', response.json())


    def test_authenticated_api_request(self):
        # Obtain access token
        token_url = reverse('token_obtain_pair')
        credentials = {
            'username': self.username,
            'password': self.password,
        }
        token_response = self.client.post(token_url, data=credentials)
        access_token = token_response.json()['access']

        # Make an authenticated API request
        api_url = reverse('api_drug')
        headers = {'HTTP_AUTHORIZATION': f'Bearer {access_token}'}
        response = self.client.get(api_url, **headers)

        self.assertEqual(response.status_code, 200)