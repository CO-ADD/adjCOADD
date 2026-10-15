#!/bin/bash
source /opt/django/etc/dj_setenv.sh

pushd /opt/django/project/${DJPROJ}

if [ "$1" = "GO" ]; then
 echo "Create Static"
#mkdir applog
#python manage.py collectstatic
fi
cd ..
if [ "$1" = "GO" ]; then
 echo "Config Static"
#mkdir static/images/mol
#chmod 775 -R static
fi

popd

