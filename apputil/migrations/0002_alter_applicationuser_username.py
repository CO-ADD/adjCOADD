# Generated by Django 3.2.6 on 2022-11-21 03:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('apputil', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='applicationuser',
            name='username',
            field=models.CharField(max_length=55, unique=True, verbose_name='user_identity_ldap'),
        ),
    ]