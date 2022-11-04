from django.dispatch import receiver
from django_auth_ldap.backend import populate_user, LDAPBackend
from .models import ApplicationUser
from aa_chem.models import Organisms, Taxonomy
from django.shortcuts import get_object_or_404
from django.db.models.signals import pre_save



@receiver(populate_user, sender=LDAPBackend)
def ldap_auth_handler(user, ldap_user, **kwargs):
    """
        Create user from ldap_user, assign to groups based the user username
    """
    user_is_staff=False
    user_is_admin=False
    user_is_appuser=True

    print("it starts")
    if ApplicationUser.objects.filter(user_id=user.username):
        print("it start....")
        appuser=ApplicationUser.objects.filter(user_id=user.username)[0]
        print((ApplicationUser.objects.filter(user_id=user.username)))
        
        try:
            if appuser.permissions=='staff':
                user_is_staff=True
                print(user_is_staff)
                user.permissions=='staff'
            elif appuser.permissions=='admin':
                user_is_staff=True
                user_is_superuser=True
            else:
                user_is_staff=False
                user_is_superuser=False           
               
        except Exception as err:
            print (err)
  
    else:
        print("Not Appuser")
        user_is_appuser=False
    
    # user.save()
    user.is_staff=user_is_staff
    user.is_admin=user_is_admin
    user.is_appuser=user_is_appuser
    user.user_id=user.username
    try:
        appuser.delete()
        print('appuser delete')
    except Exception as err:
        print(err)
    print(f'{user.username} is appuser {user.is_appuser}')
    
    user.save()
    return    








