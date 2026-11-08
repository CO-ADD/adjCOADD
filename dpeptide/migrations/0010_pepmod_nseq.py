# Manually authored migration: add nseq (number of amino acids) to
# Modification_Dictionary.

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0009_peptide_sequence_properties'),
    ]

    operations = [
        migrations.AddField(
            model_name='modification_dictionary',
            name='nseq',
            field=models.PositiveSmallIntegerField(blank=True, default=0, verbose_name='nSeq'),
        ),
    ]
