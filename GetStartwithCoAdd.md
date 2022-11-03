# Get start with Co_Add Application:

This document includes two parts: Developing and Using the Application

## Developement section

### Setting Environment

- Install conda in Linux by following the steps in [Installing on Linux] (https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).
- Install rdkit pypi overall

#### Setting test or deploy Database environment:

create a conda environment

1. create conda env <your env name> python=3.7
2. activate env
3. install the following:
   3.1 conda install -c anaconda postgresql
   3.2 conda install -c rdkit rdkit-postgresql
   3.3 pg_ctl initdb -D <path>/<your cluster name> eg. pg_ctl initdb -D ~/chem_data (meaning create cluster chem_data under user root)
   3.4 pg_ctl -D <path>/<your cluster name> start
4. using command psql postgres access to database create database and extension:
   e.g. CREATE DATABASE chembl_28;
   \q
   psql chembl_28
   create extension if not exists rdkit;

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
   pip install psycopg2
8. pip install cx-Oracle(only for oracle... need install oracle client)

9. for chemical drawing install the following:
   pip install ipython  
   pip install CairoSVG

10. django-sequences 2.7

11. pip install django-model-utils


### Module Introduction
#### Project and Applications:
project name: coadd_web contains  project setting, developing&deploying needed files and urls for all urls.
Applications:
1. aa_chem--for all chemical componds model
   contains:...
2. app-- for utility service and users model
   contains:...

 

(used techniques, connections, codes etc..)

#### Authentication System

Authentication backend authenticates against an LDAP service.
To realize this include 3 steps:

1. configuration in settings.py :
2. customize django build-in AbstractUser Model to Applicationuser Model 
3. import signals to sub-app. Building signal model to assign user with permissions(appuser, read, write, admin )
4. superusers are able to modify users and permissions

#### Data visualization

rdMolDraw2D.MolDraw2DSVG and cairosvg.svg2png




## User Guide section

### Get start

(basic functions)

#### User Data Management

(User Group and setting up permissions for superuser and staff)

#### Data Access
