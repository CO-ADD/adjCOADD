# Generated by Django 4.2.2 on 2023-07-11 09:25

from django.conf import settings
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
            name='Collab_Group',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('group_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Group ID')),
                ('group_code', models.CharField(max_length=10, unique=True, verbose_name='Group Code')),
                ('department', models.CharField(blank=True, max_length=250, verbose_name='Department')),
                ('postal_address', models.CharField(blank=True, max_length=250, verbose_name='Postal Address')),
                ('city', models.CharField(blank=True, max_length=250, verbose_name='City')),
                ('country', models.CharField(blank=True, max_length=250, verbose_name='Country')),
                ('mta_document', models.CharField(blank=True, max_length=150, verbose_name='MTA Document')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('mta_status', models.ForeignKey(blank=True, db_column='mta_status', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_MTA+', to='apputil.dictionary', verbose_name='MTA Status')),
            ],
            options={
                'db_table': 'collab_group',
            },
        ),
        migrations.CreateModel(
            name='Organisation',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('org_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='Organisation ID')),
                ('org_name', models.CharField(blank=True, max_length=250, verbose_name='Organisation')),
                ('org_code', models.CharField(blank=True, max_length=10, verbose_name='Organisation Code')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('org_type', models.ForeignKey(blank=True, db_column='org_type', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_OrgType+', to='apputil.dictionary', verbose_name='Organisation Type')),
            ],
            options={
                'db_table': 'organisation',
            },
        ),
        migrations.CreateModel(
            name='Data_Source',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('data_id', models.CharField(max_length=25, primary_key=True, serialize=False, verbose_name='Data ID')),
                ('data_name', models.CharField(blank=True, max_length=50, verbose_name='Data Name')),
                ('data_code', models.CharField(blank=True, max_length=10, verbose_name='Data Code')),
                ('description', models.CharField(blank=True, max_length=1000, verbose_name='Description')),
                ('journal', models.CharField(blank=True, max_length=50, verbose_name='Journal')),
                ('year', models.IntegerField(blank=True, verbose_name='Year')),
                ('volume', models.CharField(blank=True, max_length=50, verbose_name='Volume')),
                ('issue', models.CharField(blank=True, max_length=50, verbose_name='Issue')),
                ('page', models.CharField(blank=True, max_length=50, verbose_name='Page')),
                ('title', models.CharField(blank=True, max_length=500, verbose_name='Title')),
                ('pubmed_id', models.IntegerField(blank=True, verbose_name='PubMed ID')),
                ('doi', models.CharField(blank=True, max_length=100, verbose_name='DOI')),
                ('url', models.CharField(blank=True, max_length=100, verbose_name='URL')),
                ('authors', models.CharField(blank=True, max_length=1000, verbose_name='Authors')),
                ('patent_id', models.CharField(blank=True, max_length=500, verbose_name='Patent IDs')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('collab_group', models.ForeignKey(blank=True, db_column='collab_group', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_CollabGroup+', to='dcollab.collab_group', verbose_name='Group')),
                ('data_type', models.ForeignKey(blank=True, db_column='data_type', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_data_type', to='apputil.dictionary', verbose_name='Data Type')),
            ],
            options={
                'db_table': 'data_source',
                'ordering': ['data_name', 'data_type'],
            },
        ),
        migrations.CreateModel(
            name='Collab_User',
            fields=[
                ('astatus', models.IntegerField(db_index=True, default=0, editable=False, verbose_name='Status')),
                ('acreated_at', models.DateTimeField(editable=False, verbose_name='Created at')),
                ('aupdated_at', models.DateTimeField(editable=False, null=True, verbose_name='Updated at')),
                ('adeleted_at', models.DateTimeField(editable=False, null=True, verbose_name='Deleted at')),
                ('user_id', models.CharField(max_length=15, primary_key=True, serialize=False, verbose_name='User ID')),
                ('title', models.CharField(blank=True, max_length=15, verbose_name='Title')),
                ('first_name', models.CharField(blank=True, max_length=50, verbose_name='First Code')),
                ('last_name', models.CharField(blank=True, max_length=50, verbose_name='Last Code')),
                ('position', models.CharField(blank=True, max_length=100, verbose_name='Position')),
                ('email', models.CharField(blank=True, max_length=50, verbose_name='EMail')),
                ('phone', models.CharField(blank=True, max_length=50, verbose_name='Phone')),
                ('subscribed', models.BooleanField(blank=True, default=False, verbose_name='Newsletter')),
                ('portal_userid', models.CharField(blank=True, max_length=50, verbose_name='UserID')),
                ('portal_pw', models.CharField(blank=True, max_length=50, verbose_name='Password')),
                ('department', models.CharField(blank=True, max_length=250, verbose_name='Department')),
                ('postal_address', models.CharField(blank=True, max_length=250, verbose_name='Postal Address')),
                ('city', models.CharField(blank=True, max_length=250, verbose_name='City')),
                ('country', models.CharField(blank=True, max_length=250, verbose_name='Country')),
                ('acreated', models.ForeignKey(editable=False, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_acreated_by', to=settings.AUTH_USER_MODEL, verbose_name='Created by')),
                ('adeleted', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_adeleted_by', to=settings.AUTH_USER_MODEL, verbose_name='Deleted by')),
                ('aupdated', models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_aupdated_by', to=settings.AUTH_USER_MODEL, verbose_name='Updated by')),
                ('group', models.ForeignKey(blank=True, db_column='group', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_Group+', to='dcollab.collab_group', verbose_name='Group Membership')),
                ('organisation', models.ForeignKey(blank=True, db_column='organisation', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_Organisation+', to='dcollab.organisation', verbose_name='Organisation')),
            ],
            options={
                'db_table': 'collab_user',
            },
        ),
        migrations.AddField(
            model_name='collab_group',
            name='organisation',
            field=models.ForeignKey(blank=True, db_column='organisation', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_Organisation+', to='dcollab.organisation', verbose_name='Organisation'),
        ),
        migrations.AddField(
            model_name='collab_group',
            name='pi',
            field=models.ForeignKey(blank=True, db_column='pi', null=True, on_delete=django.db.models.deletion.DO_NOTHING, related_name='%(class)s_PI+', to='dcollab.collab_user', verbose_name='Principal Investigator'),
        ),
    ]
