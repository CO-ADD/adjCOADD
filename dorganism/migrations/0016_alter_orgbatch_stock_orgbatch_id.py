# Generated by Django 4.1.5 on 2023-04-15 12:01

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0015_alter_organism_batch_batch_id'),
    ]

    operations = [
        migrations.AlterField(
            model_name='orgbatch_stock',
            name='orgbatch_id',
            field=models.ForeignKey(db_column='orgbatch_id', editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='to_batch_organism_id', to='dorganism.organism_batch', verbose_name='OrgBatch ID'),
        ),
    ]