# Generated by Django 4.2.1 on 2023-06-30 00:56

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0011_orgbatch_image'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='orgbatch_image',
            name='orgbatch_id',
        ),
        migrations.AddField(
            model_name='orgbatch_image',
            name='img_orgbatch_id',
            field=models.ForeignKey(db_column='img_orgbatch_id', default=0, editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_img_orgbatch_id', to='dorganism.organism_batch', verbose_name='Img OrgBatch ID'),
            preserve_default=False,
        ),
    ]