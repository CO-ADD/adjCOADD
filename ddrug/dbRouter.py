from .models import Drug, VITEK_Card,VITEK_AST,VITEK_ID

Route_list=[Drug, VITEK_Card,VITEK_AST,VITEK_ID]
class DrugsRouter:
    """
    A router to control all database operations on models in the
    applications.
    """
    route_app_labels={'ddrug'}
   

    def db_for_read(self, model, **hints):
        """
        Attempts to read auth and contenttypes models go to auth_db.
        """
        if model in Route_list:
            return 'ddrug'
        return None

    def db_for_write(self, model, **hints):
        """
        Attempts to write auth and contenttypes models go to auth_db.
        """
        if model in Route_list:
            return 'ddrug'
        return None

    def allow_relation(self, obj1, obj2, **hints):
        """
        Allow relations if a model in the drugs apps is
        involved.
        """
        if (
            obj1._meta.app_label in self.route_app_labels or
            obj2._meta.app_label in self.route_app_labels
        ):
           return True
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        """
        Make sure the auth and contenttypes apps only appear in the
        'auth_db' database.
        """
        if app_label in self.route_app_labels:
            return db == 'ddrug'
        return None