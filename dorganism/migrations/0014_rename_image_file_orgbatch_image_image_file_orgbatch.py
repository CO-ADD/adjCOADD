# Generated by Django 4.2.1 on 2023-06-30 05:25

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('dorganism', '0013_alter_orgbatch_image_img_orgbatch_id'),
    ]

    operations = [
        migrations.RenameField(
            model_name='orgbatch_image',
            old_name='image_file',
            new_name='image_file_orgbatch',
        ),
    ]