from django.conf import settings

def global_settings(request):
    return {
        'DEVELOPMENT': settings.DEVELOPMENT,
        'VERSION': settings.VERSION,
        'DBHOST_NAME': settings.DBHOST_NAME,
        'DJANGO_VER': settings.DJANGO_VER,
        'PYTHON_VER': settings.PYTHON_VER,
    }