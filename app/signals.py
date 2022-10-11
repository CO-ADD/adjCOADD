from django.dispatch import receiver
from django_auth_ldap.backend import populate_user, LDAPBackend
from .models import User, Groupfilter
from django.shortcuts import get_object_or_404
from django.db.models.signals import post_save
import sys


@receiver(populate_user, sender=LDAPBackend)
def ldap_auth_handler(user, ldap_user, **kwargs):
    """
        Create user from ldap_user, assign to groups based the user email
    """
   
   
    myGroups=Groupfilter.objects.all()
    found=False
    for group in myGroups:
        if user.email in group.members: 
            user.groups.add(group) 
            print(f'add {user} in {group.name}')
            user.save(force_insert = True)
            found=True
            return 
    # if found==False:
    return None    
    #     print('user not belong to any groups')
    #     return user.delete()
from django.shortcuts import redirect
@receiver(post_save, sender=User)
def deleteuser(sender,instance,created,**kwargs):
    user=User.objects.get(username=instance.username)
    my_groups = ', '.join(map(str, user.groups.all()))
    if user.is_superuser:
        return
    else:
        try:
            if my_groups=="":
                print("no group")
                user.delete()

            #write your logic here
        #     print("User delete")
                return redirect('/accounts/login')
        # else:
        #     print(my_groups)
        except Exception as err:
            print(err)