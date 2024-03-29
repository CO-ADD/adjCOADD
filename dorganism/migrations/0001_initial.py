# Generated by Django 4.2.2 on 2023-07-11 09:21

from django.conf import settings
import django.contrib.postgres.fields
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('apputil', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Organism',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('organism_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Organism ID')),
                ('strain_ids', models.CharField(blank=True, max_length=200, verbose_name='Strain IDs')),
                ('strain_code', models.CharField(blank=True, max_length=30, verbose_name='Strain Code')),
                ('strain_panel', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=100, null=True), blank=True, null=True, size=20, verbose_name='Panel')),
                ('strain_type', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=100, null=True), blank=True, null=True, size=20, verbose_name='Type')),
                ('strain_notes', models.CharField(blank=True, max_length=1024, verbose_name='Strain Notes')),
                ('res_property', models.CharField(blank=True, max_length=1024, verbose_name='Phenotype')),
                ('gen_property', models.CharField(blank=True, max_length=1024, verbose_name='Genotype')),
                ('sero_clone', models.CharField(blank=True, max_length=126, verbose_name='MLST/Serotype')),
                ('strain_identification', models.CharField(blank=True, max_length=512, verbose_name='Strain Identification')),
                ('strain_origin', models.CharField(blank=True, max_length=512, verbose_name='Origin of Strain')),
                ('source', models.CharField(blank=True, max_length=250, verbose_name='Source')),
                ('source_code', models.CharField(blank=True, max_length=120, verbose_name='Source Code')),
                ('tax_id', models.IntegerField(default=0, verbose_name='NCBI Tax ID')),
                ('reference', models.CharField(blank=True, max_length=150, verbose_name='Reference')),
                ('mta_document', models.CharField(blank=True, max_length=150, verbose_name='MTA Document')),
                ('mta_notes', models.CharField(blank=True, max_length=512, verbose_name='MTA Notes')),
                ('collect_date', models.DateField(blank=True, null=True, verbose_name='Collection Date')),
                ('collect_region', models.CharField(blank=True, max_length=25, verbose_name='Region/City')),
                ('collect_country', models.CharField(blank=True, max_length=25, verbose_name='Country')),
                ('collect_site', models.CharField(blank=True, max_length=50, verbose_name='Site/Org')),
                ('collect_specie', models.CharField(blank=True, max_length=20, verbose_name='From Specie/Location')),
                ('collect_tissue', models.CharField(blank=True, max_length=120, verbose_name='From Tissue/Organ')),
                ('patient_diagnosis', models.CharField(blank=True, max_length=120, verbose_name='Patient Diagnosis')),
                ('patient', models.CharField(blank=True, max_length=20, verbose_name='Patient Info')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('assoc_documents', models.ManyToManyField(blank=True, db_table='org_doc', related_name='%(class)s_document', to='apputil.document', verbose_name='Douments')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('biologist', models.ForeignKey(blank=True, db_column='biologist', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_biologist', to=settings.AUTH_USER_MODEL, verbose_name='Biologist')),
                ('lab_restriction', models.ForeignKey(blank=True, db_column='lab_restriction', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_lab', to='apputil.dictionary', verbose_name='Lab')),
                ('mta_status', models.ForeignKey(blank=True, db_column='mta_status', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_mta', to='apputil.dictionary', verbose_name='MTA Status')),
            ],
            options={
                'db_table': 'organism',
                'ordering': ['organism_id'],
            },
        ),
        migrations.CreateModel(
            name='Organism_Batch',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('orgbatch_id', models.CharField(max_length=20, primary_key=True, serialize=False, verbose_name='OrgBatch ID')),
                ('batch_id', models.CharField(blank=True, max_length=12, verbose_name='Batch ID')),
                ('batch_notes', models.CharField(blank=True, max_length=500, verbose_name='Batch Notes')),
                ('qc_record', models.CharField(blank=True, max_length=150, verbose_name='QC Records')),
                ('stock_date', models.DateField(blank=True, null=True, verbose_name='Stock Date')),
                ('stock_level', models.CharField(blank=True, max_length=20, verbose_name='Stock Levels')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('biologist', models.ForeignKey(blank=True, db_column='biologist', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_biologist', to=settings.AUTH_USER_MODEL, verbose_name='Biologist')),
                ('organism_id', models.ForeignKey(db_column='organism_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_organism_id', to='dorganism.organism', verbose_name='Organism ID')),
                ('qc_status', models.ForeignKey(blank=True, db_column='qc_status', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_qc', to='apputil.dictionary', verbose_name='QC status')),
            ],
            options={
                'db_table': 'orgbatch',
                'ordering': ['orgbatch_id'],
            },
        ),
        migrations.CreateModel(
            name='Taxonomy',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('organism_name', models.CharField(max_length=100, primary_key=True, serialize=False, unique=True, verbose_name='Name')),
                ('urlname', models.SlugField(max_length=100, verbose_name='URLName')),
                ('other_names', models.CharField(blank=True, max_length=100, verbose_name='Other Names')),
                ('code', models.CharField(blank=True, max_length=15, verbose_name='Code')),
                ('tax_id', models.IntegerField(default=0, verbose_name='NCBI Tax ID')),
                ('parent_tax_id', models.IntegerField(default=0, verbose_name='NCBI Parent Tax ID')),
                ('tax_rank', models.CharField(blank=True, max_length=50, verbose_name='Taxonomy Rank')),
                ('lineage', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=60), null=True, size=30, verbose_name='Lineage')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('division', models.ForeignKey(db_column='division', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_division', to='apputil.dictionary', verbose_name='Division')),
                ('org_class', models.ForeignKey(blank=True, db_column='org_class', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_class', to='apputil.dictionary', verbose_name='Class')),
            ],
            options={
                'db_table': 'taxonomy',
                'ordering': ['organism_name'],
            },
        ),
        migrations.CreateModel(
            name='OrgBatch_Stock',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('n_created', models.IntegerField(default=0, verbose_name='#Vials created')),
                ('n_left', models.IntegerField(default=0, verbose_name='#Vials left')),
                ('stock_date', models.DateField(blank=True, null=True, verbose_name='Stock Date')),
                ('stock_note', models.CharField(blank=True, max_length=10, verbose_name='Stock Note')),
                ('location_freezer', models.CharField(blank=True, max_length=80, verbose_name='Freezer')),
                ('location_rack', models.CharField(blank=True, max_length=10, verbose_name='Rack')),
                ('location_column', models.CharField(blank=True, max_length=10, verbose_name='Column')),
                ('location_slot', models.CharField(blank=True, max_length=10, verbose_name='Slot')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('biologist', models.ForeignKey(blank=True, db_column='biologist', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_biologist', to=settings.AUTH_USER_MODEL, verbose_name='Biologist')),
                ('orgbatch_id', models.ForeignKey(db_column='orgbatch_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_orgbatch_id', to='dorganism.organism_batch', verbose_name='OrgBatch ID')),
                ('stock_type', models.ForeignKey(db_column='stock_type', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_stock', to='apputil.dictionary', verbose_name='Stock Type')),
            ],
            options={
                'db_table': 'orgbatch_stock',
                'ordering': ['orgbatch_id', 'stock_type'],
            },
        ),
        migrations.CreateModel(
            name='OrgBatch_Image',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('image_name', models.CharField(max_length=120, unique=True, verbose_name='Name')),
                ('image_file', models.ImageField(upload_to='images/orgbatch', verbose_name='Image')),
                ('image_type', models.CharField(max_length=25, verbose_name='Image Type')),
                ('image_desc', models.CharField(blank=True, max_length=140, verbose_name='Description')),
                ('image_source', models.CharField(blank=True, max_length=50, verbose_name='Source')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('orgbatch_id', models.ForeignKey(db_column='orgbatch_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_orgbatch_id', to='dorganism.organism_batch', verbose_name='OrgBatch ID')),
            ],
            options={
                'db_table': 'orgbatch_image',
                'ordering': ['orgbatch_id', 'image_name'],
            },
        ),
        migrations.CreateModel(
            name='Organism_Culture',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('media', models.CharField(blank=True, max_length=120, verbose_name='Media')),
                ('addition', models.CharField(blank=True, max_length=55, verbose_name='Addition')),
                ('atmosphere', models.CharField(blank=True, max_length=120, verbose_name='Atmosphere')),
                ('temperature', models.CharField(blank=True, max_length=25, verbose_name='Temperature')),
                ('culture_notes', models.CharField(blank=True, max_length=512, verbose_name='Notes')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('biologist', models.ForeignKey(blank=True, db_column='biologist', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_biologist', to=settings.AUTH_USER_MODEL, verbose_name='Biologist')),
                ('culture_source', models.ForeignKey(db_column='culture_source', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_culture_source', to='apputil.dictionary', verbose_name='Source')),
                ('culture_type', models.ForeignKey(db_column='culture_type', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_culture_type', to='apputil.dictionary', verbose_name='Culture Type')),
                ('organism_id', models.ForeignKey(db_column='organism_id', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_organism_id', to='dorganism.organism', verbose_name='Organism ID')),
            ],
            options={
                'db_table': 'organism_culture',
                'ordering': ['organism_id', 'culture_type', 'media'],
            },
        ),
        migrations.AddField(
            model_name='organism',
            name='organism_name',
            field=models.ForeignKey(db_column='organism_name', on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_organism_name', to='dorganism.taxonomy', verbose_name='Organism Name'),
        ),
        migrations.AddField(
            model_name='organism',
            name='oxygen_pref',
            field=models.ForeignKey(blank=True, db_column='oxygen_pref', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_oxygen', to='apputil.dictionary', verbose_name='Oxygen'),
        ),
        migrations.AddField(
            model_name='organism',
            name='pathogen_group',
            field=models.ForeignKey(blank=True, db_column='pathogen_group', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_pathogen', to='apputil.dictionary', verbose_name='Pathogen'),
        ),
        migrations.AddField(
            model_name='organism',
            name='risk_group',
            field=models.ForeignKey(blank=True, db_column='risk_group', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_risk', to='apputil.dictionary', verbose_name='Risk Group'),
        ),
        migrations.AddIndex(
            model_name='taxonomy',
            index=models.Index(fields=['org_class'], name='tax_orgclass_idx'),
        ),
        migrations.AddIndex(
            model_name='taxonomy',
            index=models.Index(fields=['tax_id'], name='tax_taxid_idx'),
        ),
        migrations.AddIndex(
            model_name='taxonomy',
            index=models.Index(fields=['division'], name='tax_div_idx'),
        ),
        migrations.AddIndex(
            model_name='taxonomy',
            index=models.Index(fields=['tax_rank'], name='tax_rnk_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_stock',
            index=models.Index(fields=['stock_type'], name='orgbstock_stype_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_stock',
            index=models.Index(fields=['location_freezer'], name='orgbstock_freezer_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_stock',
            index=models.Index(fields=['stock_date'], name='orgbstock_stdate_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_stock',
            index=models.Index(fields=['n_left'], name='orgbstock_nleft_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_image',
            index=models.Index(fields=['image_name'], name='obimg_name_idx'),
        ),
        migrations.AddIndex(
            model_name='orgbatch_image',
            index=models.Index(fields=['image_source'], name='obimg_scr_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_culture',
            index=models.Index(fields=['media'], name='orgcult_media_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_culture',
            index=models.Index(fields=['culture_type'], name='orgcult_ctype_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_batch',
            index=models.Index(fields=['organism_id', 'batch_id'], name='orgbatch_orgbatch_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_batch',
            index=models.Index(fields=['qc_status'], name='orgbatch_qc_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_batch',
            index=models.Index(fields=['stock_date'], name='orgbatch_sdate_idx'),
        ),
        migrations.AddIndex(
            model_name='organism_batch',
            index=models.Index(fields=['stock_level'], name='orgbatch_slevel_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['strain_ids'], name='org_stid_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['strain_code'], name='org_stcode_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['strain_type'], name='org_strainid_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['strain_panel'], name='org_stpanel_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['source'], name='org_source_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['tax_id'], name='org_taxid_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['risk_group'], name='org_riskgrp_idx'),
        ),
        migrations.AddIndex(
            model_name='organism',
            index=models.Index(fields=['pathogen_group'], name='org_pathgrp_idx'),
        ),
    ]
