# Generated by Django 3.2.6 on 2022-12-07 01:06

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('apputil', '0002_auto_20221205_1616'),
    ]

    operations = [
        migrations.CreateModel(
            name='AuditModel',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
        ),
        migrations.AlterField(
            model_name='applicationuser',
            name='name',
            field=models.CharField(max_length=50, primary_key=True, serialize=False, verbose_name='user name'),
        ),
        migrations.AlterField(
            model_name='applicationuser',
            name='username',
            field=models.CharField(max_length=55, unique=True, verbose_name='uq user'),
        ),
    ]