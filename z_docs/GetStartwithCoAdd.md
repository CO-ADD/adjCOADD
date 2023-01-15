# Get start with Co_Add Application:

This document includes two parts: Developing and Using the Application

## Developement section

### Setting Environment

- Install conda in Linux by following the steps in [Installing on Linux] (https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
- Install rdkit pypi overall

#### Setting test or deploy Database environment:

create a conda environment

1. create conda env <your env name>
2. activate env
3. install the following:
   3.1 conda install -c anaconda postgresql
   3.2 conda install -c rdkit rdkit-postgresql
   3.3 pg_ctl initdb -D <path>/<your cluster name> eg. pg_ctl initdb -D ~/chem_data (meaning create cluster chem_data under user root)
   3.4 pg_ctl -D <path>/<your cluster name> start
4. using command psql postgres access to database create database and extension:
   e.g. CREATE DATABASE databasename;
   \q
   psql databasename
   create extension if not exists rdkit;
   create schema...

### Setting app environment

Then create a conda environment with: `conda create -n <your-environment-name>` (this app's env is coadd_env)

#### Package installations:

Install all the following packages in the conda environment.

1. conda install -c conda-forge django
2. prerequistion for python-ldap
3. pip install python-ldap
4. pip install django-auth-ldap
5. pip install rdkit-pypi
6. pip install git+https://github.com/rdkit/django-rdkit.git
7. sudo apt-get install postgresql postgresql-contrib
   sudo apt-get install libpq-dev python3-dev
8. pip install psycopg2
9. for chemical drawing install the following:
   pip install ipython  
   pip install CairoSVG
10. django-sequences
11. pip install django-model-utils
12. pip install django-filter
13. pip install django-postgres-extra
    ...

### Module Introduction

In this section will give a structurely overview of adjCOADD firstly and then explain the main functions: authentication system and data handling(filter, crud, export and import)

#### Project and Applications:

This section provides overally introduction of functions and models of the whole project.
the adjCOADD is a Django app to screening organism database and contains 3 Applications and some util folders. A directory tree structure simply with folders is displayed below.

```
./adjCOADD
├── adjcoadd
├── apputil
│   ├── migrations
│   └── templates
│   └── apputil
├── ddrug
│   ├── migrations
│   └── templates
│   └── ddrug
│   └── drug
├── dorganism
│   ├── migrations
│   ├── templates
│   │   └── dorganism
│   │   ├── organism
│   │   │   ├── batch
│   │   │   ├── batch_stock
│   │   │   └── culture
│   │   └── taxonomy
│   └── templatetags
├── static
│   ├── css
│   ├── images
│   │   ├── app
│   │   └── brand
│   └── js
│   ├── js_utils
│   └── modal
├── templates
│   ├── registration
│   └── utils
│   └── modal
└── z_docs
```

More detailed structures will be presented in the following paragraphs and the last section.
Based on the folder structure, there are 7 main branches under project root:

1. adjcoadd -- the project core folder contains globally used constants and settings.

```
./adjCOADD
├── adjcoadd
│   ├── asgi.py
│   ├── constants.py
│   ├── __init__.py
│   ├── prod_settings.py
│   ├── routers.py
│   ├── settings.py
│   ├── urls.py
│   ├── utils_dataimport.py
│   └── wsgi.py
```

following core folder, there are 3 projects applications:

1. apputil-- prject application for utility service and users model
   models: ApplicationUser, Audit, Dictionary,
2. dorganism -- application for micro organism database
   contains:...
3. ddrug -- application for ...

The static folder is where image resources, style and interactive codes are stored under the sub-folder images, css, js respectively.

The repeated template folder followed static is contained globally used templates, e.g., base.html, nav.html...

The last folder z-doc contains help and information documents.

#### Authentication System

The authentication login system will use Ldap username and password. Once successfully log in, the applicationUser model will grant read-, write-, admin-level permissions to users.

In settings.py will set up Ldap Authentication_backend, then import Ldap and a few Ldap configuration parameters.
When user click login url, the login_user function will be called the Login_form(Django AutenticationForm) to varify if username existed in Application user datatable, then login with "django_auth_ldap.backend.LDAPBackend" backend as request user.
After login, users' permission, which has been set up previously in ApplicationUser datatable, will be granted to the login request user.

#### Data screening

For

## User Guide section

(introduction based on UI sheets)

## Full view of project - tree

```
./adjCOADD
├── adjcoadd
│   ├── asgi.py
│   ├── constants.py
│   ├── __init__.py
│   ├── prod_settings.py
│   ├── routers.py
│   ├── settings.py
│   ├── urls.py
│   ├── utils_dataimport.py
│   └── wsgi.py
├── apputil
│   ├── admin.py
│   ├── apps.py
│   ├── forms.py
│   ├── __init__.py
│   ├── migrations
│   │   ├── 0001_initial.py
│   │   ├── 0002_auto_20221205_1616.py
│   │   ├── 0003_auto_20221207_1106.py
│   │   ├── 0004_alter_dictionary_dict_desc.py
│   │   ├── 0005_alter_dictionary_acreated_alter_dictionary_adeleted_and_more.py
│   │   ├── 0006_delete_auditmodel.py
│   │   └── __init__.py
│   ├── models.py
│   ├── templates
│   │   └── apputil
│   │       ├── appUsersCreate.html
│   │       ├── appUsersDel.html
│   │       ├── appUsers.html
│   │       ├── appUsersUpdate.html
│   │       ├── appuser_tr.html
│   │       ├── dictCreate.html
│   │       ├── dictionary_del.html
│   │       ├── dictionaryUpdate.html
│   │       ├── dictList.html
│   │       └── importdata.html
│   ├── urls.py
│   ├── utils.py
│   └── views.py
├── ddrug
│   ├── admin.py
│   ├── apps.py
│   ├── dbRouter.py
│   ├── forms.py
│   ├── __init__.py
│   ├── migrations
│   │   ├── 0001_initial.py
│   │   └── __init__.py
│   ├── models.py
│   ├── templates
│   │   └── ddrug
│   │       ├── ddrug_home.html
│   │       └── drug
│   │           ├── drug_card.html
│   │           ├── drug_c.html
│   │           ├── drug_detail.html
│   │           ├── drug_d.html
│   │           ├── drug_list.html
│   │           └── drug_u.html
│   ├── tests.py
│   ├── urls.py
│   ├── utils.py
│   └── views.py
├── deploy_gunicorn.txt
├── dorganism
│   ├── admin.py
│   ├── apps.py
│   ├── dbRouter.py
│   ├── forms.py
│   ├── __init__.py
│   ├── migrations
│   │   ├── 0001_initial.py
│   │   ├── 0002_auto_20221207_0912.py
│   │   ├── 0003_taxonomy_urlname.py
│   │   ├── 0004_auto_20221207_1342.py
│   │   ├── 0005_alter_organism_options.py
│   │   ├── 0006_auto_20221215_1503.py
│   │   ├── 0007_alter_organism_options_and_more.py
│   │   └── __init__.py
│   ├── models.py
│   ├── templates
│   │   └── dorganism
│   │       ├── home.html
│   │       ├── organism
│   │       │   ├── batch
│   │       │   │   ├── batch_card.html
│   │       │   │   ├── batch_c.html
│   │       │   │   ├── batch_d.html
│   │       │   │   ├── batch_tr.html
│   │       │   │   └── batch_u.html
│   │       │   ├── batch_stock
│   │       │   │   ├── stock_c.html
│   │       │   │   ├── stock_d.html
│   │       │   │   └── stock_u.html
│   │       │   ├── culture
│   │       │   │   ├── culture_c.html
│   │       │   │   ├── culture_d.html
│   │       │   │   ├── culture_tr.html
│   │       │   │   └── culture_u.html
│   │       │   ├── organism_card.html
│   │       │   ├── organism_c.html
│   │       │   ├── organism_detail.html
│   │       │   ├── organism_d.html
│   │       │   ├── organism_list.html
│   │       │   ├── organism_u.html
│   │       │   └── organism_u_withoutname.html
│   │       └── taxonomy
│   │           ├── taxonomy_card.html
│   │           ├── taxonomy_c.html
│   │           ├── taxonomy_detail.html
│   │           ├── taxonomy_d.html
│   │           ├── taxonomy_list.html
│   │           └── taxonomy_u.html
│   ├── templatetags
│   │   ├── __init__.py
│   │   └── myapp_extras.py
│   ├── tests.py
│   ├── urls.py
│   ├── utils.py
│   └── views.py
├── .git
│   ├── branches
│   ├── config
│   ├── description
│   ├── HEAD
│   ├── hooks
│   │   ├── applypatch-msg.sample
│   │   ├── commit-msg.sample
│   │   ├── fsmonitor-watchman.sample
│   │   ├── post-update.sample
│   │   ├── pre-applypatch.sample
│   │   ├── pre-commit.sample
│   │   ├── pre-merge-commit.sample
│   │   ├── prepare-commit-msg.sample
│   │   ├── pre-push.sample
│   │   ├── pre-rebase.sample
│   │   ├── pre-receive.sample
│   │   ├── push-to-checkout.sample
│   │   └── update.sample
│   ├── index
│   ├── info
│   │   └── exclude
│   ├── logs
│   │   ├── HEAD
│   │   └── refs
│   │       ├── heads
│   │       │   └── main
│   │       └── remotes
│   │           └── origin
│   │               └── HEAD
│   ├── objects
│   │   ├── info
│   │   └── pack
│   │       ├── pack-07f70b75598eb68898d16b5b8d990a9bc5714c80.idx
│   │       └── pack-07f70b75598eb68898d16b5b8d990a9bc5714c80.pack
│   ├── packed-refs
│   └── refs
│       ├── heads
│       │   └── main
│       ├── remotes
│       │   └── origin
│       │       └── HEAD
│       └── tags
├── .gitignore
├── manage.py
├── MANIFEST.in
├── README.rst
├── requirements.txt
├── setup.cfg
├── setup.py
├── static
│   ├── css
│   │   ├── app.css
│   │   ├── bootstrap.min.css
│   │   ├── bulma.min.css
│   │   ├── custom.css
│   │   └── main.css
│   ├── images
│   │   ├── app
│   │   │   ├── avatar.png
│   │   │   └── favicon-16x16.png
│   │   └── brand
│   │       ├── CO-ADD_Baseline_White.png
│   │       ├── coaddbg01.png
│   │       ├── CO-ADD_Logo_Baseline.jpg
│   │       ├── CO-ADD logo no background.png
│   │       ├── CO-ADD_Logo_White_NoBackground_Baseline1500x375.png
│   │       ├── CO-ADD_Logo_White_NoBackground_CMYK.png
│   │       ├── CO-ADD only logo.jpg
│   │       └── logo-fullsize.png
│   └── js
│       ├── editableTablechoices.js
│       ├── editableTable.js
│       ├── importhandler.js
│       ├── index.js
│       ├── jquery-3.6.1.min.js
│       ├── js_utils
│       │   ├── ajax_submit.js
│       │   ├── ajax_update.js
│       │   ├── dataTables.js
│       │   ├── dataTables.min.js
│       │   └── getCookie.js
│       ├── modal
│       │   ├── create_batch.js
│       │   ├── create_culture.js
│       │   ├── create_drug.js
│       │   ├── create_organism.js
│       │   ├── create_stock.js
│       │   └── create_taxonomy.js
│       ├── resizebar.js
│       ├── resizebar_t.js
│       ├── search_organism_id.js
│       ├── search_organism.js
│       ├── table-editable.int.js
│       └── table-edits.min.js
├── templates
│   ├── base.html
│   ├── registration
│   │   └── login.html
│   └── utils
│       ├── alertnav.html
│       ├── card.html
│       ├── datatable_batch.html
│       ├── datatable_cultr.html
│       ├── datatable_drug.html
│       ├── datatable_general.html
│       ├── datatable_taxo.html
│       ├── editable_tr.html
│       ├── message.html
│       ├── modal
│       │   ├── create.html
│       │   ├── delete.html
│       │   └── update.html
│       ├── multichoice.html
│       ├── navbar.html
│       ├── pagination.html
│       ├── preloader.html
│       ├── search_organism.html
│       ├── search_organism_id.html
│       ├── selectAllExp.html
│       ├── sidebar.html
│       └── topbar.html
└── z_docs
    └── GetStartwithCoAdd.md

```
