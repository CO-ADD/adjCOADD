from django.apps import AppConfig


class ApputilConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'apputil'

    def ready(self):
        import apputil.signals