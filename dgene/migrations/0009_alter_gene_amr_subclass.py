# Generated by Django 4.2.2 on 2023-09-29 05:28

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dgene', '0008_alter_amr_genotype_closest_id'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='amr_subclass',
            field=models.CharField(blank=True, max_length=250, verbose_name='AMR Subclass'),
        ),
    ]
