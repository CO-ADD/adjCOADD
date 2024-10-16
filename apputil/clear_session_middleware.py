from django.urls import resolve
from django.urls.exceptions import Resolver404

class ClearSessionMiddleware:
    # pass
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        # Get the current and last visited views
       # If resolver_match exists, get the current_view
     
        sessionkeyslist=list(request.session.keys())
        try:
            current_view = resolve(request.path_info).url_name
        except Resolver404:
            current_view = None
        last_view = request.session.get('last_view')
            # If the current view is different from the last visited view, and both are in clear_session_views, clear the session data
        if last_view and last_view != current_view and current_view!='dataexport' and current_view!='pivoted-table':
            for i in list(sessionkeyslist):
                if 'cached_queryset' in i:              
                    del request.session[i]
 
            # Update the last visited view
        request.session['last_view'] = current_view
        # Process the request
        response = self.get_response(request)
                
        return response
