from django.dispatch import receiver
from django_auth_ldap.backend import populate_user, LDAPBackend
from .models import ApplicationUser
from django.shortcuts import get_object_or_404
from django.db.models.signals import pre_save

@receiver(populate_user, sender=LDAPBackend)
def ldap_auth_handler(user, ldap_user, **kwargs):
    """
        Create user from ldap_user, assign to groups based the user username
    """

    user.save()
    print("it starts")
    appuser= ApplicationUser.objects.filter(user_id=user.username)
    print (appuser)
    if appuser:
        try:
            if appuser.permissions=='staff':
                user.is_staff=True
            elif appuser.permissions=='admin':
                user.is_staff=True
                user.is_admin=True
            else:
                user.is_staff=False
                user.is_admin=False
               
        except Exception as err:
            print (err)
  
    else:
        user.is_appuser=False

    return    



# def pre_save(sender, instance, created, **kwargs):
    