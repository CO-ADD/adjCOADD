# Generated by Django 4.2.2 on 2023-09-29 03:06

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('dgene', '0002_alter_gene_amr_subclass'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='source',
            field=models.CharField(blank=True, max_length=150, verbose_name='Source'),
        ),
    ]