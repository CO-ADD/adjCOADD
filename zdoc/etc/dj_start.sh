#!/bin/bash -
source /opt/django/etc/dj_setenv.sh

gunicorn -c /opt/django/etc/dj_gunicorn.py --daemon 

