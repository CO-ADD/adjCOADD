from django.db import models

# Create your models here.
from django.db import models
from django.contrib.auth.models import AbstractUser, Group, Permission
from django.contrib.postgres.fields import ArrayField
from django.urls import reverse

# Create your models here.



class User(AbstractUser):
   
    uq_id = models.IntegerField(null=True)
    organisation = models.CharField(max_length=250, null=True)
    department = models.CharField(max_length=250, null=True)
    avatar = models.ImageField(null=True, upload_to='media')
    createddate=models.DateTimeField(null=True,verbose_name = "Created on", auto_now_add=True, editable=False)

    def get_short_name(self):
        # Returns the short name for the user.
        fStr = self.first_name
        try:
            iStr =fStr[0]
            return f"{iStr}.{self.last_name}"
        except Exception as err:
            print (err)
            return

    # def get_username(self):
    #     # Returns the short name for the user.
    #     return f"{self.id})"

    def __str__(self):
        return f"{self.username} ({self.email})"

    def is_admin(user):
        return user.groups.filter(name='admin').exists()
   
class Groupfilter(Group):
    filtername= models.CharField(max_length=250, null=True, unique=True)
    description=models.CharField(max_length=250, null=True)
    members=ArrayField(
        models.CharField(max_length=100, blank=True, null=True), size=20, null=True, blank=True
    )
    
    def get_absolute_url(self):
        return reverse('usermanage')

    class Meta:
        permissions = (('importdata', 'change data'),)
  
