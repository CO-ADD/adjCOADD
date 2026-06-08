import sys, os
#import types
import datetime
import numpy as np
import pandas as pd
import re

import psycopg2
from applib.external.sql_connector import SqlConnector
from sshtunnel import SSHTunnelForwarder

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
# PostgreSQL
#
# For sshtunnel - install paramiko=3.5.1 
#                 paramiko=4.x does not work with sshtunnel 0.4
#
#-----------------------------------------------------------------------------
class PostgreSQL(SqlConnector):

    def open(self, username, password, hostname, database, 
             ssh_user=None, ssh_password=None, ssh_private_key=None,  verbose=1):


        self.user = username
        self.password = password
        self.hostname = hostname
        self.port = 5432
        self.database = database
        
        try:
            self.config = {
                    'user': self.user,
                    'password': self.password,
                    'host': self.hostname,
                    'port' : self.port,
                    'dbname': self.database
                    #'raise_on_warnings': True
                    }
            
            # Direct Connetion    
            if ssh_user is None:
                self.tunnel = None

                #self.db = psycopg2.connect(host=hostname, dbname=database, user=username, password=password)
                self.db = psycopg2.connect(h**self.config)
                if self.verbose>0:
                    logger.info(f"[PostgreSQL] {self.database}@{self.hostname}")
            
            # SSH-Tunnel Connetion    
            elif ssh_user is not None and ssh_password is not None:

                print(f" {ssh_user} {ssh_password} {ssh_private_key}")
                self.tunnel = SSHTunnelForwarder(
                    hostname,
                    ssh_username=ssh_user,
                    ssh_password=ssh_password,
                    #ssh_pkey=ssh_private_key,
                    remote_bind_address=('127.0.0.1', 5432),
                    local_bind_address=('localhost',6543), 
                )
                # Start the tunnel
                self.tunnel.start()
                print(f" SSH Tunnel started ")
                
                self.config['host'] = self.tunnel.local_bind_host
                self.config['port'] = self.tunnel.local_bind_port
                self.db = psycopg2.connect(**self.config)
                if self.verbose>0:
                    logger.info(f"[PostgreSQL] {self.database}@{self.hostname} [SSH]")
                
        except psycopg2.OperationalError as e:
            logger.error(f"[PostgreSQL] Connection Error - {e}")
            raise

        self.verbose=verbose
        self.sql_type = "PostgreSQL"
        
        self.cursor = self.db.cursor()
        #self.cursor.execute('SET GLOBAL max_allowed_packet=500*1024*1024')
        
    def exec(self, sql, bindvars=None, commit=False):
        try:
            self.cursor.execute(sql, bindvars)
        except psycopg2.OperationalError as e:
            logger.error(f"[PostgreSQL] Execute Error {e}")
            raise
        if commit:
            self.db.commit()

    def execmany(self, sql, bindvars=None, commit=False):
        try:
            self.cursor.executemany(sql, bindvars)
        except psycopg2.OperationalError as e:
            logger.error(f"[PostgreSQL] Execute Many Error {e}")
            raise
        if commit:
            self.db.commit()

    def close(self):
        try:
            #self.cursor.close()
            self.db.close()
            if self.tunnel is not None:
                self.tunnel.close()
                
        except psycopg2.connector.Error as e:
            pass
        if self.verbose>0:
            logger.info(f"[PostgreSQL] {self.database}@{self.hostname} [Closed] ")

    def sqlbindvar(self,str):
        return("%("+str+")s")

    def gen_sqltyping(self,var,varname=None):

        if isinstance(var,(int,float)):
            if varname is None:
                rSql = "to_number("+ str(var) + ")"
            else:
                #rSql = "to_number("+ self.sqlbindvar(varname) + ")"
                rSql = self.sqlbindvar(varname)
        elif isinstance(var,str):
            if varname is None:
                rSql = self.sqlquote(var)
            else:
                rSql = self.sqlbindvar(varname)
        elif isinstance(var,datetime.datetime):
            if varname is None:
                rSql = "to_date("+ var.strftime("%Y-%m-%d %H-%M-%S") + ",'DD-MM-YYYY HH-MI-SS')"
            else:
                # rSql = "to_date("+ self.sqlbindvar(varname) + ",'DD-MM-YYYY HH-MI-SS')"
                rSql = self.sqlbindvar(varname)
        elif isinstance(var,datetime.date):
            if varname is None:
                rSql = "to_date('"+ var.strftime("%Y-%m-%d") + "','YYYY-MM-DD')"
            else:
                #rSql = "to_date("+ self.sqlbindvar(varname) + ",'DD-MM-YYYY')"
                rSql = self.sqlbindvar(varname)
        else:
            if varname is None:
                rSql = var
            else:
                rSql = self.sqlbindvar(varname)
        return(rSql)
