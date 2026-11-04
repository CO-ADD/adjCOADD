# Manually authored migration: rename bilm -> biln and widen sequence
# notation fields (seq/biln/helm) from CharField(200) to TextField.

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0003_remove_aminoacid_dictionary_aa_id_idx_and_more'),
    ]

    operations = [
        # Drop indexes on the renamed bilm fields (they become TextField, which
        # cannot carry a models.Index).
        migrations.RemoveIndex(
            model_name='peptide',
            name='pep_bilm_idx',
        ),
        migrations.RemoveIndex(
            model_name='peptide_isoform',
            name='pepiso_bilm_idx',
        ),

        # Rename bilm -> biln on all three models.
        migrations.RenameField(
            model_name='peptide',
            old_name='bilm',
            new_name='biln',
        ),
        migrations.RenameField(
            model_name='peptide_isoform',
            old_name='bilm',
            new_name='biln',
        ),
        migrations.RenameField(
            model_name='modification_dictionary',
            old_name='bilm',
            new_name='biln',
        ),

        # Widen the sequence notation fields to TextField.
        migrations.AlterField(
            model_name='peptide',
            name='seq',
            field=models.TextField(blank=True, verbose_name='Seq'),
        ),
        migrations.AlterField(
            model_name='peptide',
            name='biln',
            field=models.TextField(blank=True, verbose_name='BILN'),
        ),
        migrations.AlterField(
            model_name='peptide',
            name='helm',
            field=models.TextField(blank=True, verbose_name='HELM'),
        ),
        migrations.AlterField(
            model_name='peptide_isoform',
            name='seq',
            field=models.TextField(blank=True, verbose_name='Seq'),
        ),
        migrations.AlterField(
            model_name='peptide_isoform',
            name='biln',
            field=models.TextField(blank=True, verbose_name='BILN'),
        ),
        migrations.AlterField(
            model_name='peptide_isoform',
            name='helm',
            field=models.TextField(blank=True, verbose_name='HELM'),
        ),
        migrations.AlterField(
            model_name='modification_dictionary',
            name='seq',
            field=models.TextField(blank=True, verbose_name='Seq'),
        ),
        migrations.AlterField(
            model_name='modification_dictionary',
            name='biln',
            field=models.TextField(blank=True, verbose_name='BILN'),
        ),
        migrations.AlterField(
            model_name='modification_dictionary',
            name='helm',
            field=models.TextField(blank=True, verbose_name='HELM'),
        ),
    ]
