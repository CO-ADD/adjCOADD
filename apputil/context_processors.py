from django.conf import settings

def global_settings(request):
    return {
        'DEVELOPMENT': settings.DEVELOPMENT,
        'VERSION': settings.VERSION,
        'HOST_NAME': settings.HOST_NAME,
    }