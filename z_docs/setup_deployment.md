#
# Deployment on Production using nginx and gunicron
#
# gunicorn - python WSGI HTTP server
# nginx - Advanced load balancer, Web server and Reserve Proxy
#

#-------------------------------------------------------------------------------------
# Install gunicorn ino conda env
#-------------------------------------------------------------------------------------
dj41py310> conda install gunicorn  #gunicorn-20.1.0  

#-------------------------------------------------------------------------------------
# Install nginx
#-------------------------------------------------------------------------------------
root> yum install epel-release
root> yum install nginx  #nginx 1:1.20.1-10.el7

#-------------------------------------------------------------------------------------
# start nginx service
#-------------------------------------------------------------------------------------
root> systemctl start nginx
root> systemctl status nginx
 
[root@imb-coadd-work ~]# systemctl status nginx.service
        nginx.service - The nginx HTTP and reverse proxy server
        Loaded: loaded (/usr/lib/systemd/system/nginx.service; disabled; vendor preset: disabled)
        Process: 5377 ExecStartPre=/usr/bin/rm -f /run/nginx.pid (code=exited, status=0/SUCCESS)
        Process: 5385 ExecStartPre=/usr/sbin/nginx -t (code=exited, status=0/SUCCESS)
        Process: 5388 ExecStart=/usr/sbin/nginx (code=exited, status=0/SUCCESS)
        Active: active (running) since Wed 2022-12-21 14:05:17 AEST; 2s ago
        Main PID: 5390 (nginx)
            Tasks: 3
        CGroup: /system.slice/nginx.service
                ├─5390 nginx: master process /usr/sbin/nginx
                ├─5391 nginx: worker process
                └─5392 nginx: worker process

        Dec 21 14:05:17 imb-coadd-work systemd[1]: Starting The nginx HTTP and reverse proxy server...
        Dec 21 14:05:17 imb-coadd-work nginx[5385]: nginx: the configuration file /etc/nginx/nginx.conf syntax is ok
        Dec 21 14:05:17 imb-coadd-work nginx[5385]: nginx: configuration file /etc/nginx/nginx.conf test is successful
        Dec 21 14:05:17 imb-coadd-work systemd[1]: Started The nginx HTTP and reverse proxy server.


root> systemctl eneable nginx  #to enable start nginx on boot-up
root> systemctl restart nginx
root> systemctl stop nginx

#-------------------------------------------------------------------------------------
# configuration of nginx service
#-------------------------------------------------------------------------------------
# /etc/nginx/
root> more /etc/nginx/nginx.conf
        # For more information on configuration, see:
        #   * Official English Documentation: http://nginx.org/en/docs/
        #   * Official Russian Documentation: http://nginx.org/ru/docs/

        user nginx;
        worker_processes auto;
        error_log /var/log/nginx/nginx_error.log;
        pid /run/nginx.pid;

        # Load dynamic modules. See /usr/share/doc/nginx/README.dynamic.
        include /usr/share/nginx/modules/*.conf;

        events {
            worker_connections 1024;
        }

        http {
            log_format  main  '$remote_addr - $remote_user [$time_local] "$request" '
                            '$status $body_bytes_sent "$http_referer" '
                            '"$http_user_agent" "$http_x_forwarded_for"';

            access_log  /var/log/nginx/access.log  main;

            sendfile            on;
            tcp_nopush          on;
            tcp_nodelay         on;
            keepalive_timeout   65;
            types_hash_max_size 4096;

            include             /etc/nginx/mime.types;
            default_type        application/octet-stream;

            # Load modular configuration files from the /etc/nginx/conf.d directory.
            # See http://nginx.org/en/docs/ngx_core_module.html#include
            # for more information.
            include /etc/nginx/conf.d/*.conf;

            server {
->              listen       8080;
->              listen       [::]:8080;
                server_name  _;
                root         /usr/share/nginx/html;

                # Load configuration files for the default server block.
->              include /etc/nginx/default.d/*.conf;

                error_page 404 /404.html;
                location = /404.html {
                }

                error_page 500 502 503 504 /50x.html;
                location = /50x.html {
                }
            }

        # Settings for a TLS enabled server.
        #
        #    server {
        #        listen       443 ssl http2;
        #        listen       [::]:443 ssl http2;
        #        server_name  _;
        #        root         /usr/share/nginx/html;
        #
        #        ssl_certificate "/etc/pki/nginx/server.crt";
        #        ssl_certificate_key "/etc/pki/nginx/private/server.key";
        #        ssl_session_cache shared:SSL:1m;
        #        ssl_session_timeout  10m;
        #        ssl_ciphers HIGH:!aNULL:!MD5;
        #        ssl_prefer_server_ciphers on;
        #
        #        # Load configuration files for the default server block.
        #        include /etc/nginx/default.d/*.conf;
        #
        #        error_page 404 /404.html;
        #            location = /40x.html {
        #        }
        #
        #        error_page 500 502 503 504 /50x.html;
        #            location = /50x.html {
        #        }
        #    }

        }

#-------------------------------------------------------------------------------------
# Allow port access for services
#-------------------------------------------------------------------------------------
semanage port -a -t http_port_t -p tcp 8090

#-------------------------------------------------------------------------------------
# /etc/nginx/conf.d/django.conf
#-------------------------------------------------------------------------------------

        #
        # Gunicorn configuration for Django
        #

        server {
            server_name         imb-coadd-work.uq.edu.au;
            listen              8008;
            listen              [::]:8008;

            location / {
            ## semanage fcontext -a -t httpd_t "/run/gunicorn(/.*)?"
            ## restorecon -r /run
->          #proxy_pass       http://unix:/opt/temp/gunicorn.sock;  # did not work due to selinux permision issues 
-?          proxy_pass        http://localhost:8009;
            proxy_set_header  Host $host;
            proxy_set_header  X-Forwarded-For $proxy_add_x_forwarded_for;
            }

->          location /static {
->          autoindex         on;
->          alias             /opt/django/project/static;
            }

->          location /uploads {
->          autoindex         on;
->          alias             /opt/django/var/uploads;
            }
        }


#-------------------------------------------------------------------------------------
<h2> Configuration of gunicorn </h2>
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
# /opt/django/etc/dj_gunicorn.py
#-------------------------------------------------------------------------------------
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
        workers = 4

        # The socket to bind
->      #bind = 'unix:/opt/django/var/gunicorn.sock' # did not work due to selinux permision issues 
->      #bind = 'unix:/opt/temp/gunicorn.sock'       # did not work due to selinux permision issues 
->      bind = "0.0.0.0:8009"

        # Restart workers when code changes (development only!)
        reload = True

        # Write access and error info to /var/log
        logfile   = "/opt/django/var/log/gunicorn.log"
        accesslog = "/opt/django/var/log/gunicorn_access.log"
        errorlog  = "/opt/django/var/log/gunicorn_error.log"
        worker_tmp_dir = "/opt/django/var/tmp"

        # Redirect stdout/stderr to log file
        capture_output = True

        # PID file so you can easily fetch process ID
        pidfile = "/opt/django/var/gunicorn.pid"

        # Daemonize the Gunicorn process (detach & enter background)
        daemon = True

        #user="nginx"
        #group="nginx"

        timeout=120


#-------------------------------------------------------------------------------------
<h2> Configuration of Firewall </h2>
#-------------------------------------------------------------------------------------

# /usr/lib/firewalld/services/

#-------------------------------------------------------------------------------------
django.xml
#-------------------------------------------------------------------------------------
<?xml version="1.0" encoding="utf-8"?>
<service>
  <short>Django</short>
  <description>Django User ports</description>
  <port protocol="tcp" port="8008"/>
  <port protocol="tcp" port="5000"/>
</service>

#-------------------------------------------------------------------------------------
postgresql.xml
#-------------------------------------------------------------------------------------
<?xml version="1.0" encoding="utf-8"?>
<service>
  <short>PostgreSQL</short>
  <description>PostgreSQL Database Server</description>
  <port protocol="tcp" port="5432"/>
</service>


root> firewall-cmd –reload

root> firewall-cmd --get-active-zones
	public
	  interfaces: ens192 

root> firewall-cmd --zone=public --add-service=<newservices>
root> firewall-cmd --zone=public --permanent --add-service=<newservices>
root> firewall-cmd --zone=public --list-services


#-------------------------------------------------------------------------------------
<h2> Configuration of Systemctl </h2>
#-------------------------------------------------------------------------------------

# /usr/lib/systemd/system/

#-------------------------------------------------------------------------------------
django.service
#-------------------------------------------------------------------------------------
[Unit]
Description=Django WebApplication Server
After=network.target

[Service]
Type=forking

User=django
Group=django

ExecStart=/opt/django/etc/dj_start.sh
ExecStop=/opt/django/etc/dj_stop.sh

StandardOutput=null

# Give a reasonable amount of time for the server to start up/shut down
TimeoutSec=300

#-------------------------------------------------------------------------------------
postgresql.service
#-------------------------------------------------------------------------------------
[Unit]
Description=PostgreSQL database server
After=network.target

[Service]
Type=forking

User=postgres
Group=postgres

# Disable OOM kill on the postmaster
OOMScoreAdjust=-1000

ExecStart=/opt/postgres/etc/pg_start.sh
ExecStop=/opt/postgres/etc/pg_stop.sh
ExecReload=/opt/postgres/etc/pg_reload.sh

# Give a reasonable amount of time for the server to start up/shut down
TimeoutSec=300




root> systemctl daemon-reload
root> systemctl list-unit-files --type service | grep <service>
root> systemctl enable <service>

root> systemctl status <service>
