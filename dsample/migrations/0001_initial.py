# Generated by Django 4.2.13 on 2024-07-02 20:50

from django.conf import settings
import django.contrib.postgres.fields
import django.core.validators
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('apputil', '0004_alter_applicationuser_options'),
        ('dcollab', '0001_initial'),
        ('dchem', '0001_initial'),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('sample_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Sample ID')),
                ('batch_id', models.CharField(blank=True, default='00', max_length=12, validators=[django.core.validators.RegexValidator('^[0-9a-zA-Z]*$', 'Only alphanumeric characters are allowed.')], verbose_name='Batch ID')),
                ('batch_notes', models.CharField(blank=True, max_length=500, verbose_name='Batch Notes')),
                ('sample_code', models.CharField(blank=True, max_length=150, verbose_name='Sample Code')),
                ('sample_name', models.CharField(blank=True, max_length=250, verbose_name='Sample Name')),
                ('sample_desc', models.CharField(blank=True, max_length=512, verbose_name='Sample Description')),
                ('previous_ids', models.CharField(blank=True, max_length=100, verbose_name='Previous IDs')),
                ('salt_code', models.CharField(blank=True, max_length=120, verbose_name='Salts')),
                ('full_mw', models.FloatField(blank=True, default=0, verbose_name='Full MW')),
                ('full_mf', models.CharField(blank=True, max_length=100, verbose_name='Full MF')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('sample_type', models.ForeignKey(blank=True, db_column='sample_type', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_sample_type', to='apputil.dictionary', verbose_name='Sample Type')),
                ('structure_id', models.ForeignKey(blank=True, db_column='structure_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_structure_id', to='dchem.chem_structure', verbose_name='Structure ID')),
            ],
            options={
                'db_table': 'sample',
                'ordering': ['sample_id'],
            },
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('project_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Project ID')),
                ('old_project_id', models.CharField(max_length=15, unique=True, verbose_name='Old Project ID')),
                ('project_name', models.CharField(blank=True, max_length=50, verbose_name='Project Name')),
                ('owner_user', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=25, null=True), blank=True, null=True, size=10, verbose_name='User')),
                ('source', models.CharField(blank=True, max_length=250, verbose_name='Source')),
                ('source_code', models.CharField(blank=True, max_length=120, verbose_name='Source Code')),
                ('reference', models.CharField(blank=True, max_length=150, verbose_name='Reference')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('owner_group', models.ForeignKey(blank=True, db_column='owner_group', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_ownergroup', to='dcollab.collab_group', verbose_name='Group')),
                ('project_class', models.ForeignKey(blank=True, db_column='project_class', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_projectclass', to='apputil.dictionary', verbose_name='Class')),
            ],
            options={
                'db_table': 'project',
                'ordering': ['project_id'],
            },
        ),
        migrations.CreateModel(
            name='Convert_ProjectID',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('ora_project_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Old Project ID')),
                ('project_id', models.CharField(max_length=15, verbose_name='Project ID')),
                ('project_name', models.CharField(blank=True, max_length=50, verbose_name='Project Name')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
            ],
            options={
                'db_table': 'convert_projectid',
                'ordering': ['ora_project_id'],
            },
        ),
        migrations.CreateModel(
            name='Convert_CompoundID',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('ora_compound_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Old Compound ID')),
                ('compound_id', models.CharField(max_length=15, verbose_name='Compound ID')),
                ('compound_code', models.CharField(blank=True, max_length=120, verbose_name='Code')),
                ('compound_name', models.CharField(blank=True, max_length=120, verbose_name='Name')),
                ('project_id', models.CharField(max_length=15, verbose_name='Project ID')),
                ('sample_type', models.CharField(max_length=10, verbose_name='Project ID')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
            ],
            options={
                'db_table': 'convert_compoundid',
                'ordering': ['ora_compound_id'],
            },
        ),
        migrations.CreateModel(
            name='COADD_Compound',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('compound_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Compound ID')),
                ('compound_code', models.CharField(blank=True, max_length=120, verbose_name='Code')),
                ('compound_name', models.CharField(blank=True, max_length=120, verbose_name='Name')),
                ('compound_desc', models.CharField(blank=True, max_length=150, verbose_name='Comment')),
                ('old_compound_id', models.CharField(max_length=15, unique=True, verbose_name='Old ID')),
                ('reg_smiles', models.CharField(blank=True, max_length=2048, verbose_name='Reg Smiles')),
                ('reg_mw', models.FloatField(blank=True, default=0, verbose_name='Reg MW')),
                ('reg_mf', models.CharField(blank=True, max_length=100, verbose_name='Reg MF')),
                ('reg_structure', models.CharField(blank=True, max_length=2048, verbose_name='Reg Structure')),
                ('reg_amount', models.FloatField(blank=True, default=0, verbose_name='Reg Amount')),
                ('reg_amount_unit', models.CharField(blank=True, max_length=100, verbose_name='Reg Amount Unit')),
                ('reg_volume', models.FloatField(blank=True, default=0, verbose_name='Reg Volume')),
                ('reg_volume_unit', models.CharField(blank=True, max_length=100, verbose_name='Reg Volume Unit')),
                ('reg_conc', models.FloatField(blank=True, default=0, verbose_name='Reg Conc')),
                ('reg_conc_unit', models.CharField(blank=True, max_length=100, verbose_name='Reg Conc Unit')),
                ('reg_solvent', models.FloatField(blank=True, default=0, verbose_name='Reg Solvent')),
                ('prep_date', models.DateField(blank=True, null=True, verbose_name='Prepared')),
                ('stock_amount', models.FloatField(blank=True, default=0, verbose_name='Stock Amount')),
                ('stock_amount_unit', models.CharField(blank=True, max_length=100, verbose_name='Stock Amount Unit')),
                ('std_status', models.CharField(blank=True, max_length=10, verbose_name='Std Status')),
                ('std_action', models.CharField(blank=True, max_length=120, verbose_name='Std Action')),
                ('std_process', models.CharField(blank=True, max_length=120, verbose_name='Std Process')),
                ('std_smiles', models.CharField(blank=True, max_length=2048, verbose_name='Std Smiles')),
                ('std_mw', models.FloatField(blank=True, default=0, verbose_name='Std MW')),
                ('std_mf', models.CharField(blank=True, max_length=100, verbose_name='Std MF')),
                ('cpoz_sn', models.CharField(blank=True, max_length=25, verbose_name='CpOz SN')),
                ('cpoz_id', models.CharField(blank=True, max_length=25, verbose_name='CpOz Lib ID')),
                ('coadd_id', models.CharField(blank=True, max_length=25, verbose_name='CO-ADD ID')),
                ('chembl_id', models.CharField(blank=True, max_length=25, verbose_name='ChEMBL ID')),
                ('spark_id', models.CharField(blank=True, max_length=25, verbose_name='SPARK ID')),
                ('pub_status', models.CharField(blank=True, max_length=10, verbose_name='Pub Status')),
                ('pub_date', models.DateField(blank=True, editable=False, null=True, verbose_name='Published')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('compound_type', models.ForeignKey(blank=True, db_column='compound_type', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_compound_type', to='apputil.dictionary', verbose_name='Type')),
                ('project_id', models.ForeignKey(blank=True, db_column='project_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_project_id', to='dsample.project', verbose_name='Project ID')),
                ('sample_id', models.ForeignKey(blank=True, db_column='sample_id', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_sample_id', to='dsample.sample', verbose_name='Sample ID')),
            ],
            options={
                'db_table': 'coadd_sample',
                'ordering': ['sample_id'],
            },
        ),
        migrations.AddIndex(
            model_name='sample',
            index=models.Index(fields=['sample_name'], name='sample_name_idx'),
        ),
        migrations.AddIndex(
            model_name='sample',
            index=models.Index(fields=['sample_code'], name='sample_code_idx'),
        ),
        migrations.AddIndex(
            model_name='sample',
            index=models.Index(fields=['sample_type'], name='sample_type_idx'),
        ),
        migrations.AddIndex(
            model_name='sample',
            index=models.Index(fields=['full_mw'], name='sample_fmw_idx'),
        ),
        migrations.AddIndex(
            model_name='sample',
            index=models.Index(fields=['salt_code'], name='sample_salt_idx'),
        ),
        migrations.AddIndex(
            model_name='project',
            index=models.Index(fields=['project_name'], name='prj_pname_idx'),
        ),
        migrations.AddIndex(
            model_name='project',
            index=models.Index(fields=['old_project_id'], name='prj_opid_idx'),
        ),
        migrations.AddIndex(
            model_name='convert_projectid',
            index=models.Index(fields=['project_id'], name='wprj_pid_idx'),
        ),
        migrations.AddIndex(
            model_name='convert_projectid',
            index=models.Index(fields=['project_name'], name='wprj_pname_idx'),
        ),
        migrations.AddIndex(
            model_name='convert_compoundid',
            index=models.Index(fields=['compound_id'], name='wcmpd_cid_idx'),
        ),
        migrations.AddIndex(
            model_name='convert_compoundid',
            index=models.Index(fields=['compound_code'], name='wcmpd_ccode_idx'),
        ),
        migrations.AddIndex(
            model_name='convert_compoundid',
            index=models.Index(fields=['sample_type'], name='wcmpd_stype_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['compound_name'], name='coadd_name_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['compound_code'], name='coadd_code_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['compound_type'], name='coadd_type_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['old_compound_id'], name='coadd_ocid_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['std_status'], name='coadd_sst_idx'),
        ),
        migrations.AddIndex(
            model_name='coadd_compound',
            index=models.Index(fields=['pub_status'], name='coadd_pst_idx'),
        ),
    ]