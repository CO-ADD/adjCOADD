# Manually authored migration: add computed peptide properties (pI, GRAVY,
# charge at pH 7) to Peptide and Peptide_Isoform.

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0008_alter_modification_dictionary_fasta_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='peptide',
            name='isoelectric_point',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=6, null=True, verbose_name='pI'),
        ),
        migrations.AddField(
            model_name='peptide',
            name='gravy',
            field=models.DecimalField(blank=True, decimal_places=3, max_digits=8, null=True, verbose_name='Hydrophobicity (GRAVY)'),
        ),
        migrations.AddField(
            model_name='peptide',
            name='charge_at_ph7',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=6, null=True, verbose_name='Charge @ pH 7'),
        ),
        migrations.AddField(
            model_name='peptide_isoform',
            name='isoelectric_point',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=6, null=True, verbose_name='pI'),
        ),
        migrations.AddField(
            model_name='peptide_isoform',
            name='gravy',
            field=models.DecimalField(blank=True, decimal_places=3, max_digits=8, null=True, verbose_name='Hydrophobicity (GRAVY)'),
        ),
        migrations.AddField(
            model_name='peptide_isoform',
            name='charge_at_ph7',
            field=models.DecimalField(blank=True, decimal_places=2, max_digits=6, null=True, verbose_name='Charge @ pH 7'),
        ),
    ]
