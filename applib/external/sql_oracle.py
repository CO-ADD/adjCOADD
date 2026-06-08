import sys, os
#import types
import datetime
import numpy as np
import pandas as pd
import re

#import cx_Oracle as ora
import oracledb as ora

from applib.external.sql_connector import SqlConnector 

import logging
logger = logging.getLogger(__name__)


def IntConv(value):
    return int(value)  # or whatever is needed to convert from numpy.int64 to an integer

def InputTypeHandler(cursor, value, num_elements):
    if isinstance(value, np.int64):
        return cursor.var(int, arraysize=num_elements, inconverter=IntConv)

#-----------------------------------------------------------------------------
# Oracle
#-----------------------------------------------------------------------------
class Oracle(SqlConnector):

    def open(self, username, password, hostname, port, servicename,verbose=1):
        self.verbose = verbose
        try:
            self.dsn = ora.makedsn(hostname,port,service_name=servicename)
            self.db = ora.connect(user=username,password=password, dsn=self.dsn)
        except ora.DatabaseError as e:
            err = e.args
            logger.error(f"[Oracle] Connection Error {err}")
            raise
        self.sql_type = "Oracle"
        self.cursor = self.db.cursor()
        self.user = username
        self.password = password
        self.hostname = hostname
        self.port = port
        self.servicename = servicename
        # self.cursor.execute('select max(count(*))  from v$open_cursor group by sid')
        self.cursor.execute("Select sys_context('USERENV','SID') From Dual")
        self.sid = self.cursor.fetchone()[0]
        if self.verbose>0:
            logger.info("[Oracle] {:s}@{:s} [SID: {:s}]".format(self.user,self.servicename,self.sid))
        self.cursor.inputtypehandler = InputTypeHandler

    def close(self):
        try:
            self.cursor.close()
            self.db.close()
        except ora.DatabaseError as e:
            pass
        if self.verbose>0:
            logger.info("[Oracle] {:s}@{:s} [Closed] ".format(self.user,self.servicename))

    def exec(self, sql, bindvars=None, commit=False):
        try:
            if (bindvars == None):
                self.cursor.execute(sql)
            else:
                xbindvars = self.fix_sqlvars(sql,bindvars)
                self.cursor.execute(sql,xbindvars)
        except ora.DatabaseError as e:
            print(sql)
            if bindvars != None:
                print(bindvars)
            #logger.error("[Oracle] Database Error {:s} {:s}".format(errObj.code,errObj.message))
            #logger.error("[Oracle] Database Error {:s}".format(sql))
            #if bindvars != None:
            #    logger.error("[Oracle] Database Error {:s}".format(bindvars))
            raise
        except ora.IntegrityError as e:
            errObj = e.args
            logger.error("[Oracle] Execute Error {:s} {:s}".format(errObj.code,errObj.message))
            raise
        if commit:
            self.db.commit()

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
