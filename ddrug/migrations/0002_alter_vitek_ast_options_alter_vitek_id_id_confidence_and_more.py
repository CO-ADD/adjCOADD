# Generated by Django 4.1.1 on 2023-01-23 10:22

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ddrug', '0001_initial'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='vitek_ast',
            options={'ordering': ['drug_id', 'bp_source', 'card_barcode']},
        ),
        migrations.AlterField(
            model_name='vitek_id',
            name='id_confidence',
            field=models.CharField(blank=True, max_length=120, verbose_name='ID Confidence'),
        ),
        migrations.AlterField(
            model_name='vitek_id',
            name='id_probability',
            field=models.CharField(blank=True, max_length=120, verbose_name='ID Probability'),
        ),
        migrations.AddIndex(
            model_name='vitek_ast',
            index=models.Index(fields=['bp_source'], name='vast_bpsource_idx'),
        ),
        migrations.AddIndex(
            model_name='vitek_id',
            index=models.Index(fields=['card_barcode'], name='vid_barcode_idx'),
        ),
    ]