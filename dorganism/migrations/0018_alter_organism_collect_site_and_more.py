# Generated by Django 4.2.1 on 2023-07-07 07:02

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0017_alter_organism_batch_qc_status_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='organism',
            name='collect_site',
            field=models.CharField(blank=True, max_length=50, verbose_name='Site/Org'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='collect_tissue',
            field=models.CharField(blank=True, max_length=120, verbose_name='From Tissue/Organ'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='mta_notes',
            field=models.CharField(blank=True, max_length=512, verbose_name='MTA Notes'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='patient_diagnosis',
            field=models.CharField(blank=True, max_length=120, verbose_name='Patient Diagnosis'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='sero_clone',
            field=models.CharField(blank=True, max_length=126, verbose_name='MLST/Serotype'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='strain_identification',
            field=models.CharField(blank=True, max_length=512, verbose_name='Strain Identification'),
        ),
        migrations.AlterField(
            model_name='organism',
            name='strain_origin',
            field=models.CharField(blank=True, max_length=512, verbose_name='Origin of Strain'),
        ),
    ]