from django.dispatch import receiver
from django_auth_ldap.backend import populate_user, LDAPBackend
from .models import User, Groupfilter
from django.shortcuts import get_object_or_404

@receiver(populate_user, sender=LDAPBackend)
def ldap_auth_handler(user, ldap_user, **kwargs):
    """
        Create user from ldap_user, assign to groups based the user username
    """
   
   
    user.save()
    myGroups=Groupfilter.objects.all()
    found=False
    for group in myGroups:
        if user.username in group.members: 
            user.groups.add(group) 
            print(f'add {user} in {group.name}')
            found=True
            return 
    
    return    
