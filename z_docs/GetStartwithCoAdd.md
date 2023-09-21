# Get start with CO-ADD Application:

### Production Conda
- Download Conda from https://docs.conda.io/en/main/miniconda.html -> Miniconda3-latest-Linux-x86_64.sh
- Install conda as django user 
   django> sh Miniconda3-latest-Linux-x86_64.sh 
      /opt/django/conda

### Production [lnx] CentOS 7.- Setting Environment - dj42py311

- (base) django> conda update -n base -c defaults conda

- (base) django> conda create -n dj42py311 django django-rdkit <br>
   creates base enviroment with django 4.2.2, python 3.11.3

- (base) django> conda activate <b>dj42py311</b>

- (dj42py311) django> <b>conda</b> install ..<br>
      django-filter django-pandas psycopg2 python-ldap requests djangorestframework <br>
      pdfplumber python-magic openpyxl  ipython CairoSVG<br>
      rdkit -c rdkit<br>
      gunicorn<br>
      postgresql<br>

   [win] django-filter django-pandas psycopg2 requests djangorestframework <br>
   [win] openpyxl ipython CairoSVG<br>
   [win] rdkit -c rdkit  <br>
    conda install libpython m2w64-toolchain <br>
   [win for import data] cx_oracle mysql-connector-python scipy [for zLib]<br>

- (dj42py311) django> <b>pip</b> install ..<br>
   django-model-utils django-sequences django-formtools<br>
   django-auth-ldap <br>
   djangorestframework-simplejwt<br>
   clamd<br>
   django-dbbackup<br>
   django-crontab <br>

   [win] django-model-utils django-sequences django-formtools
   [win] djangorestframework-simplejwt
   [win] python-magic pdfplumber

   [win] get python_ldap-3.4.3-cp311-cp311-win_amd64.whl from https://github.com/cgohlke/python-ldap-build/releases<br>
   [win] pip install ldap-3.4.3-cp311-cp311-win_amd64.whl<br>
   [win] django-auth-ldap <br>


   [win for import data] synergy [for zLib]
### Development Local [win] Win10 .- Setting Environment - dj42py311

- (base) > conda update -n base -c defaults conda

- (base) > conda create -n dj42py311 django django-rdkit<br>
   creates base enviroment with django 4.2.2, python 3.11.3

- (base) > conda activate <b>dj42py311</b>

- (dj42py311) > <b>conda</b> install ..<br>
   django-filter django-pandas psycopg2 requests djangorestframework <br>
   openpyxl ipython CairoSVG<br>
   rdkit -c rdkit  <br>
   (libpython m2w64-toolchain) - not sure if needed 


- (dj42py311) > <b>pip</b> install ..<br>
   django-model-utils django-sequences django-formtools<br>
   djangorestframework-simplejwt<br>
   python-magic pdfplumber

- (dj42py311) > <b>pip</b> install django-auth-ldap<br>
   installed VC Build Tools from https://visualstudio.microsoft.com/visual-cpp-build-tools/ <br>
      Visual C++ build tools <br>
      To avoid: Pip error: Microsoft Visual C++ 14.0 is required

   from https://stackoverflow.com/questions/15070063/installing-python-ldap-in-a-virtualenv-on-windows/58235938#58235938<br>
   get python_ldap-3.4.3-cp311-cp311-win_amd64.whl from https://github.com/cgohlke/python-ldap-build/releases<br>
      for py311 on bit (AMD64)<br>

      (dj42py311) > python -c 'import sys; print(sys.version)'<br>
      3.11.3 | packaged by conda-forge | (main, Apr  6 2023, 08:50:54) [MSC v.1934 64 bit (AMD64)]

   pip install ldap-3.4.3-cp311-cp311-win_amd64.whl

   pip install django-auth-ldap 

   or try https://github.com/etianen/django-python3-ldap

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

##### Way 1

- step1. run: conda env create -f environment.yml
- step2. install python-ldap (https://www.python-ldap.org/en/python-ldap-3.4.3/installing.html)
- step3. pip install django-auth-ldap
- step4. pip install git+https://github.com/rdkit/django-rdkit.git
- step5. (Linux)sudo apt-get install libmagic1; (windows) download .dll and .mgc 3 files from https://github.com/pidydx/libmagicwin64/tree/master, copy to the project root.
- step6. pip install jinjia2

##### Way 2 :

Manually install all the following packages in the conda environment.

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
13. pip install pdfplumber == 0.7.6
14. pip install
15. pip install djangorestframework
16. pip install djangorestframework_simplejwt
17. pip install requests
    (this is just for test purpose)
18. pip install openpyxl
19. pip install django-formtools
20. pip install clamd
21. sudo apt-get install libmagic1
22. pip install python-magic
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

1. apputil-- project application for utility service and users model
   models: ApplicationUser, Audit, Dictionary,
2. dorganism -- application for micro organism database
   contains:...
3. ddrug -- application for ...

The static folder is where image resources, style and interactive codes are stored under the sub-folder images, css, js respectively.

The repeated template folder followed static is contained globally used templates, e.g., base.html, nav.html...

The last folder z-doc contains help and information documents.

#### Authentication System

The authentication login system will use Ldap username and password. Once successfully log in, the applicationUser model will be given read-, write-, admin-level permissions to users.

In settings.py will set up Ldap Authentication_backend, then import Ldap and a few Ldap configuration parameters.
When user click login url, the login_user function will be called the Login_form(Django AutenticationForm) to varify if username existed in Application user datatable, then login with "django_auth_ldap.backend.LDAPBackend" backend as request user.
After login, users' permission, which has been set up previously in ApplicationUser datatable, will be granted to the login request user.

#### Data screening

All data tables are sorting table type. And in project settings.py, the Django-Filter packages has been implemented and it is used to get filtered queryset from database for the main models(Organism, Taxonomy, Drug, Vitek_Card).
To display one to many relationship, the child table will be displayed under the parent table's detail view.

#### Data CRUD

To create data in database, A bootstrap modal form containing "django ModelForm (for create)" is used to create a new entry for each model.
For browsing all the data in main tables, there are list- and cardview for each model displaying entries per page via django pagination.
To update entries editable tables and modal form are both used. Models updated by editable tables are: Organism, Batch, Culture, Applicationusers. Models updated by modal forms are: Stock, Taxonomy.

## User Guide section

### Login

The adjCOADD is started with login using ldap authentication.

<img src="https://github.com/CO-ADD/adjCOADD/blob/main/static/images/app/CoAdd_Login.png" />

### Homepage

Once login successfully, will direct to the homepage:

<img src="https://github.com/CO-ADD/adjCOADD/blob/main/static/images/app/CoAdd_Home.png" />

With main datatables' access.

### Files upload
