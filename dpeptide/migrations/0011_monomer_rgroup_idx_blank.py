# Manually authored migration: allow blank (None) items inside the
# rgroup_idx / attachment_idx ArrayFields.

import django.contrib.postgres.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0010_pepmod_nseq'),
    ]

    operations = [
        migrations.AlterField(
            model_name='monomer',
            name='rgroup_idx',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(blank=True, null=True), blank=True, null=True, size=4, verbose_name='R-Group Idx'),
        ),
        migrations.AlterField(
            model_name='monomer',
            name='attachment_idx',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.IntegerField(blank=True, null=True), blank=True, null=True, size=4, verbose_name='Attachment Idx'),
        ),
    ]
