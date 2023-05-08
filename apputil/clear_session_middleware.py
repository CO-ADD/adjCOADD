from django.urls import reverse

class ClearSessionMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        
        # List the views for which you want to clear the session
        clear_session_views = [
            reverse('org_card'),
            reverse('org_list'),
            reverse('taxo_list'),
            reverse('taxo_card'),
            reverse('taxo_card'),
            reverse('wgs-fastqc_list'),
            reverse('wgs-checkm_list'),
        ]

        # Clear the session if the user is navigating to any of the listed views
        if request.path in clear_session_views:
            if 'cached_queryset' in request.session:
                del request.session['cached_queryset']
                
        return response
