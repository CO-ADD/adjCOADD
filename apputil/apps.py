from django.apps import AppConfig


class ApptilConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'apputil'

    def ready(self):
        pass