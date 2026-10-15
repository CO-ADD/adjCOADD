#
# Gunicorn settings for Django Production on imb-coadd
#
import multiprocessing

# Django WSGI application path in pattern MODULE_NAME:VARIABLE_NAME
chdir = "/opt/django/project/adjCOADD"
wsgi_app = "adjcoadd.wsgi:application"

# The granularity of Error log outputs
#loglevel = 'info'
loglevel = "debug"

# The number of worker processes for handling requests
# workers = multiprocessing.cpu.count()*2+1
workers = 9

# The socket to bind
#bind = 'unix:/opt/django/var/gunicorn.sock'
#bind = 'unix:/opt/temp/gunicorn.sock'
#bind = "0.0.0.0:8009"
bind = "127.0.0.1:8009"

# Restart workers when code changes (development only!)
reload = True

# Write access and error info to /var/log
logfile   = "/opt/django/var/log/gunicorn.log"
accesslog = "/opt/django/var/log/gunicorn_access.log"
errorlog  = "/opt/django/var/log/gunicorn_error.log"
worker_tmp_dir = "/opt/django/var/tmp"

# Redirect stdout/stderr to log file
#capture_output = True

# PID file so you can easily fetch process ID
pidfile = "/opt/django/var/gunicorn.pid"

# Daemonize the Gunicorn process (detach & enter background)
#daemon = True

#user="nginx"
#group="nginx"

timeout=120

