# Generated by Django 4.1.3 on 2023-05-15 06:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0002_alter_organism_gen_property_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='organism',
            name='gen_property',
            field=models.CharField(blank=True, max_length=1024, verbose_name='Genetic Property'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='res_property',
            field=models.CharField(blank=True, max_length=1024, verbose_name='Susceptibility '),
        ),
        migrations.AlterField(
            model_name='organism',
            name='strain_notes',
            field=models.CharField(blank=True, max_length=1024, verbose_name='Strain Notes'),
        ),
    ]