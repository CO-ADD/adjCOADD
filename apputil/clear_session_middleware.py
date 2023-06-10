from django.urls import reverse

class ClearSessionMiddleware:
    # pass
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        response = self.get_response(request)
        # Get the current and last visited views
        current_view = request.resolver_match.url_name
        last_view = request.session.get('last_view')
        # List the views 
        clear_session_views = [
            reverse('userslist'),
            reverse('dict_view'),

            reverse('api_drug'),
            reverse('api_vitekast'),
            reverse('drug_list'),
            reverse('drug_card'),
            reverse('vitekcard_list'),
            reverse('vitekast_list'),
            reverse('vitekid_list'),
            reverse('mic_coadd_list'),
            reverse('mic_coadd_card'),
            reverse('mic_pub_list'),
            reverse('mic_pub_card'),
            reverse('breakpoint_list'),

            reverse('gene_list'),
            reverse('id_pub_list'),
            reverse('sequence_list'),
            reverse('wgs-fastqc_list'),
            reverse('wgs-checkm_list'),

            reverse('org_card'),
            reverse('org_list'),
            reverse('taxo_list'),
            reverse('taxo_card'),
        ]

        # If the current view is different from the last visited view, and both are in clear_session_views, clear the session data
        if last_view and last_view != current_view and current_view in clear_session_views and last_view in clear_session_views:
            if 'cached_queryset' in request.session:
                del request.session['cached_queryset']
        
        # Update the last visited view
        request.session['last_view'] = current_view
                
        return response
