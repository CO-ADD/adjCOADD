"""
Django settings for adjCOADD project.

Generated by 'django-admin startproject' using Django 4.1.1.

For more information on this file, see
https://docs.djangoproject.com/en/4.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/4.1/ref/settings/
"""
import os
from pathlib import Path
# import cx_Oracle
# cx_Oracle.init_oracle_client(lib_dir="/opt/oracle/instantclient_21_7")

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent

MEDIA_ROOT=os.path.join(BASE_DIR, 'uploads')
MEDIA_URL=('uploads/')


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/4.1/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'django-insecure-_fzrv(t#j+r4y)7s$nm=v!qt=+!@vs(2-=z)ls(h^$ozyj!$g^'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ["*"]


# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',   
    'django_rdkit',
    'django_filters',
    "sequences.apps.SequencesConfig",
    "django.contrib.postgres",
    "psqlextra",
    'apputil.apps.ApptilConfig',
    'dorganism.apps.DorganismConfig',
    'ddrug.apps.DdrugConfig',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'adjcoadd.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [
            BASE_DIR/'templates',
        #    BASE_DIR/,
            ],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'adjcoadd.wsgi.application'


# Database
# https://docs.djangoproject.com/en/4.1/ref/settings/#databases

DATABASES = {
    'default': {
        # 'ENGINE': 'django.db.backends.postgresql_psycopg2',
        
        "ENGINE": "psqlextra.backend",
        'OPTIONS':{'options': '-c search_path=apputil,public'},
        'NAME': 'orgdb',
        'USER': 'orgdb', #os.environ.get('db_user'),
        'PASSWORD':'orgdb',
        'HOST': 'imb-coadd-work.imb.uq.edu.au',
        'PORT': '5432',
    },
    'dorganism': {
        "ENGINE": "psqlextra.backend",
        'OPTIONS':{'options': '-c search_path=dorganism,apputil'},
        'NAME': 'orgdb',
        'USER': 'orgdb', #os.environ.get('db_user'),
        'PASSWORD': 'orgdb',
        'HOST': 'imb-coadd-work.imb.uq.edu.au',
        'PORT': '5432',
    },

    'ddrug': {
        "ENGINE": "psqlextra.backend",
        'OPTIONS':{'options': '-c search_path=ddrug,dorganism,apputil,public'},
        'NAME': 'orgdb',
        'USER': 'orgdb', #os.environ.get('db_user'),
        'PASSWORD': 'orgdb',
        'HOST': 'imb-coadd-work.imb.uq.edu.au',
        'PORT': '5432',
    }


    # 'default': {
    #     "ENGINE": "psqlextra.backend",
    #     'OPTIONS':{'options': '-c search_path=apputil,public'},
    #     'NAME': 'orgdb_local',
    #     # 'USER': 'postgres', #os.environ.get('db_user'),
    #     # 'PASSWORD':os.environ.get('db_password', 'password'),
       
    #     'HOST': 'Localhost',
    #     'PORT': '5432',
    # },
    # 'dorganism': {
    #     "ENGINE": "psqlextra.backend",
    #     'OPTIONS':{'options': '-c search_path=dorganism,apputil'},
    #     'NAME': 'orgdb_local',
    #     'TEST': {
    #         'NAME': 'mytestdatabase2',
    #     },
    #     # 'USER': 'postgres', #os.environ.get('db_user'),
    #     # 'PASSWORD': os.environ.get('db_password','password'),
    #     'HOST': 'Localhost',
    #     'PORT': '5432',
    # },
    #   'ddrug': {
    #     "ENGINE": "psqlextra.backend",
    #     'OPTIONS':{'options': '-c search_path=ddrug,dorganism,apputil'},
    #     'NAME': 'orgdb_local',
        
    #     # 'USER': 'postgres', #os.environ.get('db_user'),
    #     # 'PASSWORD': os.environ.get('db_password','password'),
    #     'HOST': 'Localhost',
    #     'PORT': '5432',
    # }

   


}

DATABASE_ROUTERS = ['dorganism.dbRouter.DrugsRouter','ddrug.dbRouter.DrugsRouter']
# Password validation
# https://docs.djangoproject.com/en/4.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/4.1/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/4.1/howto/static-files/

STATIC_URL = 'static/'
STATICFILES_DIRS=[BASE_DIR/"static",]
# Default primary key field type
# https://docs.djangoproject.com/en/4.1/ref/settings/#default-auto-field

AUTH_USER_MODEL = 'apputil.ApplicationUser'
DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

import ldap
from django_auth_ldap.config import LDAPSearch, GroupOfNamesType, LDAPGroupQuery, PosixGroupType


# LOGIN_REDIRECT_URL='/'
# LOGOUT_REDIRECT_URL='/accounts/login'
#LDAP AUthen
AUTHENTICATION_BACKENDS = [
    "django_auth_ldap.backend.LDAPBackend", 
    "django.contrib.auth.backends.ModelBackend",
    ]
AUTH_LDAP_SERVER_URI = "ldap://ldap.uq.edu.au"

AUTH_LDAP_BIND_DN = ""

AUTH_LDAP_BIND_PASSWORD = ""

AUTH_LDAP_USER_SEARCH = LDAPSearch("ou=people,o=The University of Queensland,c=au", ldap.SCOPE_SUBTREE, "(uid=%(user)s)")

# AUTH_LDAP_GROUP_SEARCH = LDAPSearch("ou=groups,o=The University of Queensland,c=au", ldap.SCOPE_SUBTREE, "(objectClass=groupOfNames)")

# AUTH_LDAP_GROUP_TYPE = GroupOfNamesType()

 

# AUTH_LDAP_CONNECTION_OPTIONS = {

#     ldap.OPT_DEBUG_LEVEL: 0,

#     ldap.OPT_REFERRALS: 0,

# }

# AUTH_LDAP_ALWAYS_UPDATE_USER=False
 

# AUTH_LDAP_USER_ATTR_MAP = {

#     "first_name": "givenName",

#     "last_name": "sn",

#     "email": "mail",

#     "username": "uid",

#     "uq_id": "prism",

#     "department": "ou",

# }


# ===============================================================

Strain_Type_choices=(("a", "testa"),("b","testb"),("c", "testc"))#querysetToChoiseList_Dictionaries(Dictionaries, Organisms.Choice_Dictionaries['Strain_Type'])

# Security Setting
CSRF_COOKIE_SECURE=True