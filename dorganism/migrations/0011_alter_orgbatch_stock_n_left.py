# Generated by Django 4.1.5 on 2023-03-23 06:33

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0010_alter_orgbatch_stock_n_left'),
    ]

    operations = [
        migrations.AlterField(
            model_name='orgbatch_stock',
            name='n_left',
            field=models.IntegerField(default=0, verbose_name='#Vials left'),
        ),
    ]