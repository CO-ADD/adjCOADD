# Generated by Django 4.2.1 on 2023-07-06 12:03

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0016_alter_orgbatch_image_options_and_more'),
        ('apputil', '0005_alter_document_doc_desc_alter_document_doc_source_and_more'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='image',
            name='acreated',
        ),
        migrations.RemoveField(
            model_name='image',
            name='adeleted',
        ),
        migrations.RemoveField(
            model_name='image',
            name='aupdated',
        ),
        migrations.RemoveField(
            model_name='document',
            name='id',
        ),
        migrations.AddField(
            model_name='dictionary',
            name='dict_app',
            field=models.CharField(default=1, max_length=30, verbose_name='Application'),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='document',
            name='doc_id',
            field=models.CharField(default=1, max_length=15, primary_key=True, serialize=False, verbose_name='Doc ID'),
            preserve_default=False,
        ),
        migrations.AlterField(
            model_name='document',
            name='doc_name',
            field=models.CharField(max_length=120, verbose_name='Name'),
        ),
        migrations.AddIndex(
            model_name='applicationlog',
            index=models.Index(fields=['log_status'], name='log_status_idx'),
        ),
        migrations.DeleteModel(
            name='Image',
        ),
    ]