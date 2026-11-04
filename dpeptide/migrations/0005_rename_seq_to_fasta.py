# Manually authored migration: rename seq -> fasta to reflect that the field
# stores a FASTA one-letter amino-acid sequence.

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0004_rename_bilm_to_biln_and_make_sequence_fields_text'),
    ]

    operations = [
        migrations.RenameField(
            model_name='peptide',
            old_name='seq',
            new_name='fasta',
        ),
        migrations.RenameField(
            model_name='peptide_isoform',
            old_name='seq',
            new_name='fasta',
        ),
        migrations.RenameField(
            model_name='modification_dictionary',
            old_name='seq',
            new_name='fasta',
        ),
    ]
