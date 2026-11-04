# Manually authored migration: Monomer library model.

import django.contrib.postgres.fields
import django.contrib.postgres.indexes
import django.db.models.deletion
import django_rdkit.models.fields
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0005_rename_seq_to_fasta'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Monomer',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('monomer_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Monomer ID')),
                ('symbol', models.CharField(max_length=30, unique=True, verbose_name='Symbol')),
                ('monomer_name', models.CharField(blank=True, max_length=150, verbose_name='Name')),
                ('monomer_abbr', models.CharField(blank=True, max_length=30, verbose_name='Abbreviation')),
                ('monomer_type', models.CharField(blank=True, max_length=15, verbose_name='Type')),
                ('monomer_subtype', models.CharField(blank=True, max_length=20, verbose_name='SubType')),
                ('nat_analog', models.CharField(blank=True, max_length=10, verbose_name='Natural Analog')),
                ('pdb_name', models.CharField(blank=True, max_length=10, verbose_name='PDB Name')),
                ('mol', django_rdkit.models.fields.MolField(blank=True, null=True, verbose_name='Structure')),
                ('molfile', models.TextField(blank=True, verbose_name='Mol Block')),
                ('rgroups', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=10, null=True), blank=True, null=True, size=4, verbose_name='R-Groups')),
                ('rgroup_idx', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), blank=True, null=True, size=4, verbose_name='R-Group Idx')),
                ('attachment_idx', django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(null=True), blank=True, null=True, size=4, verbose_name='Attachment Idx')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
            ],
            options={
                'db_table': 'monomer',
                'ordering': ['symbol'],
            },
        ),
        migrations.AddIndex(
            model_name='monomer',
            index=models.Index(fields=['symbol'], name='monomer_symbol_idx'),
        ),
        migrations.AddIndex(
            model_name='monomer',
            index=models.Index(fields=['monomer_type'], name='monomer_type_idx'),
        ),
        migrations.AddIndex(
            model_name='monomer',
            index=django.contrib.postgres.indexes.GistIndex(fields=['mol'], name='monomer_mol_idx'),
        ),
    ]
