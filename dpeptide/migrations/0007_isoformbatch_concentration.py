# Manually authored migration: add concentration/amount fields to Isoform_Batch.

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dpeptide', '0006_monomer'),
    ]

    operations = [
        migrations.AddField(
            model_name='isoform_batch',
            name='conc',
            field=models.DecimalField(blank=True, decimal_places=3, default=0, max_digits=12, verbose_name='Concentration'),
        ),
        migrations.AddField(
            model_name='isoform_batch',
            name='conc_unit',
            field=models.ForeignKey(blank=True, db_column='conc_unit', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_conc_unit', to='apputil.dictionary', verbose_name='Concentration Unit'),
        ),
        migrations.AddField(
            model_name='isoform_batch',
            name='init_amount',
            field=models.DecimalField(decimal_places=3, default=0, max_digits=12, verbose_name='Init Amount'),
        ),
        migrations.AddField(
            model_name='isoform_batch',
            name='init_amount_unit',
            field=models.ForeignKey(blank=True, db_column='init_amount_unit', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_init_amount_unit', to='apputil.dictionary', verbose_name='Init Amount Unit'),
        ),
    ]
