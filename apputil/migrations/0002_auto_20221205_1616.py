# Generated by Django 3.2.6 on 2022-12-05 06:16

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('apputil', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='applicationuser',
            name='name',
            field=models.CharField(max_length=50, primary_key=True, serialize=False, verbose_name='user'),
        ),
        migrations.AlterField(
            model_name='applicationuser',
            name='username',
            field=models.CharField(max_length=55, unique=True, verbose_name='uquser'),
        ),
    ]
