from django.dispatch import receiver
from django_auth_ldap.backend import populate_user, LDAPBackend
from .models import User, Groupfilter, ApplicationUser
from django.shortcuts import get_object_or_404

@receiver(populate_user, sender=LDAPBackend)
def ldap_auth_handler(user, ldap_user, **kwargs):
    """
        Create user from ldap_user, assign to groups based the user username
    """
   
   
    user.save()
    print("it starts")
    appuser= ApplicationUser.objects.filter(user_id=user.username)
    print (appuser)
    if appuser!=None:
        try:
            if appuser.permissions=='staff':
                user.role='staff'
                user.is_staff=True
            elif appuser.permissions=='admin':
                
                user.role='admin'
                user.is_admin=True
            else:
                user.is_staff=False
                user.is_admin=False
                user.role='general'
        except Exception as err:
            print (err)
    # for group in myGroups:
    #     if user.username in group.members: 
    #         user.groups.add(group) 
    #         print(f'add {user} in {group.name}')
    #         found=True
    #         return 
    else:
        user.role="None"

    return    
