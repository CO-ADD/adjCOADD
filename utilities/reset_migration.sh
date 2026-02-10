#!/bin/bash

APP_LIST01=("dcollab" "dchem" "dorganism" "dcell" "dpeptide" "ddrug" "dgene" "dsample" ) 
APP_LIST02=("dscreen" "dplate" "dsummary")

if [[ "${1}" == "RESET" ]]
then
  echo = 111 ===================
  for app in ${APP_LIST01[@]}; do	  
    echo [${app}] migrations
    #rm -rf ${app}/migrations
  done

  for app in ${APP_LIST02[@]}; do
    echo [${app}] migrations
    #rm -rf ${app}/migrations
  done

  echo = 222 ===================
  for app in ${APP_LIST01[@]}; do
    echo [${app}] Makemigrations	  
    python manage.py makemigrations ${app}
    echo ------------------
  done
  echo = 333 ===================
  for app in ${APP_LIST01[@]}; do
    echo [${app}] Fake Migrite 
    python manage.py migrate ${app} --database ${app} --fake
    echo ------------------
  done

  for app in ${APP_LIST02[@]}; do
    echo [${app}] Makemigrations
    #python manage.py makemigrations ${app}
    echo ------------------
  done
fi

if [[ "${1}" == "TEST" ]]
then
  for app in ${APP_LIST02[@]}; do
    echo [${app}] migrations
    #rm -rf ${app}/migrations
  done

  for app in ${APP_LIST02[@]}; do
    echo [${app}] Fake Migrite
    #python manage.py makemigrations ${app}
    #python manage.py migrate ${app} --database ${app} --fake
    echo ------------------
  done
fi
