import sys, os
#import types
import datetime
import numpy as np
import pandas as pd
import re

import mysql.connector
from applib.external.sql_connector import SqlConnector 

import logging
logger = logging.getLogger(__name__)

#-----------------------------------------------------------------------------
# MySQL
#-----------------------------------------------------------------------------
class MySQL(SqlConnector):

    def open(self, username, password, hostname, database, verbose = 1):
        try:
            self.config = {
                    'user': username,
                    'password': password,
                    'host': hostname,
                    'database': database,
                    'port' : 3306,
                    'raise_on_warnings': True
                    }

            self.db = mysql.connector.connect(**self.config)

        except mysql.connector.Error as e:
            logger.error("[MySQL] Connection Error - {}".format(e))
            raise
        self.verbose = verbose
        self.sql_type = "MySQL"
        self.cursor = self.db.cursor(buffered=True)
        #self.cursor.execute('SET GLOBAL max_allowed_packet=500*1024*1024')
        self.user = username
        self.password = password
        self.hostname = hostname
        self.port = 3306
        self.database = database
        self.SID = self.db.connection_id
        if self.verbose>0:
            logger.info("[MySQL] {:s}@{:s} [SID: {:d}] ".format(self.database,self.hostname,self.SID))

    def exec(self, sql, bindvars=None, commit=False):
        try:
            self.cursor.execute(sql, bindvars)
        except mysql.connector.Error as e:
            logger.error("[MySQL] Execute Error {}".format(e))
            raise
        if commit:
            self.db.commit()

    def execmany(self, sql, bindvars=None, commit=False):
        try:
            self.cursor.executemany(sql, bindvars)
        except mysql.connector.Error as e:
            logger.error("[MySQL] Execute Many Error {}".format(e))
            raise
        if commit:
            self.db.commit()

    def close(self):
        try:
            self.cursor.close()
            self.db.close()
        except mysql.connector.Error as e:
            pass
        if self.verbose>0:
            logger.info("[MySQL] {:s}@{:s} [Closed] ".format(self.database,self.hostname))

    def sqlbindvar(self,str):
        return(":"+str)

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
                # rSql = "to_date(:"+ varname + ",'DD-MM-YYYY HH-MI-SS')"
                rSql = self.sqlbindvar(varname)
        elif isinstance(var,datetime.date):
            if varname is None:
                rSql = "to_date('"+ var.strftime("%Y-%m-%d") + "','YYYY-MM-DD')"
            else:
                #rSql = "to_date(:"+ varname + ",'DD-MM-YYYY')"
                rSql = self.sqlbindvar(varname)
        else:
            if varname is None:
                rSql = var
            else:
                rSql = self.sqlbindvar(varname)
        return(rSql)

