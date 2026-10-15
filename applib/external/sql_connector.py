"""
SQL Connector Class with Oracle and MySQL options
"""
#-----------------------------------------------------------------------------
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
#-----------------------------------------------------------------------------
import sys, os
#import types
import datetime
import numpy as np
import pandas as pd
import re

import logging
logger = logging.getLogger(__name__)
__version__ = "1.1"
__author__ = "J.Zuegg UQ/CO-ADD"

#-----------------------------------------------------------------------------
# Dict utilities
#-----------------------------------------------------------------------------
def conv_dict_upperkeys(d):
    n = {}
    for k,v in d.items():
        n[k.upper()] = v
    return(n)
        #d[k.upper()] = d.pop(k)

def merge_dict(a,b, overwrite=False, keeponly=False):
    if overwrite:
        n = {**a, **b}
    else:
        n = {**b, **a}
    if keeponly:
        sc = list(set(a).intersection(b))
        return({ k: n[k] for k in sc })
    else:
        return(n)
    return

def remove_None_dict(dict):
    dictf = {k: v for k, v in dict.items() if v is not None}
    #dictf = {k: v for k, v in dictf.items() if not np.isnan(v)}
    dict.clear()
    dict.update(dictf)

def remove_List_dict(dict):
    dictf = {k: v for k, v in dict.items() if not isinstance(v,list)}
    dict.clear()
    dict.update(dictf)

#-----------------------------------------------------------------------------
# SQL Connector Class
#-----------------------------------------------------------------------------
class SqlConnector(object):

    def __init__(self):
        self.db = None
        self.cursor = None
        self.sql_type = None
        self.verbose = 1

    def open(self, username, password, hostname, verbose):
        raise NotImplementedError("Each SqlConnector is responsible for its own open method.")

    def exec(self, sql, bindvars=None, commit=False):
        raise NotImplementedError("Each SqlConnector is responsible for its own exec method.")

    def execmany(self, sql, bindvars=None, commit=False):
        raise NotImplementedError("Each SqlConnector is responsible for its own execmany method.")

    def close(self):
        raise NotImplementedError("Each SqlConnector is responsible for its own close method.")

    def nCount(self, sql, bindvars=None):
        """
        General method to return the first single column in a Select statement
          - to be used with "Select count(1) From ... "
        """
        self.exec(sql,bindvars=bindvars)
        new_Count = self.cursor.fetchone()[0]
        return new_Count

    def nget_single_dict(self, sql, bindvars=None, upCase=True):
        self.exec(sql,bindvars=bindvars)
        sqlDict = {}
        columns = [i[0] for i in self.cursor.description]
        if upCase:
            sqlColumns = [s.upper() for s in columns]
        else:
            sqlColumns = columns

        row = self.cursor.fetchone()
        if row:
            for col in sqlColumns:
                sqlDict[col] = row[sqlColumns.index(col)]

        return(sqlDict)

    def get_single_dict(self, sql, bindvars=None, upCase=True):
        self.exec(sql,bindvars=bindvars)
#        if self.cursor.rowcount > 0: DOES NOT Work for Oracle until cursor.fetchone() or cursor.fetchall()
        columns = [i[0] for i in self.cursor.description]
        new_dict = dict()
        for row in self.cursor:
            for col in columns:
                if upCase:
                    new_dict[col.upper()] = row[columns.index(col)]
                else:
                    new_dict[col] = row[columns.index(col)]
            break
        return new_dict
    
    def nget_dict_list(self, sql, bindvars=None, upCase=True, batchsize=100):
  
        self.exec(sql,bindvars=bindvars)
        sqlLst = []
        columns = [i[0] for i in self.cursor.description]
        if upCase:
            sqlColumns = [s.upper() for s in columns]
        else:
            sqlColumns = columns

        while True:
            rows = self.cursor.fetchmany(batchsize)
            if not rows:
                break
            for row in rows:
                rowDict = {}
                for col in sqlColumns:
                    rowDict[col] = row[sqlColumns.index(col)]
                sqlLst.append(rowDict)
        return(sqlLst)
            
    def get_dict_list(self, sql, bindvars=None,upCase=False,lowCase=True):
        """
        Create a list, each item contains a dictionary outlined like:
        { "col1_name" : col1_data }
        Each item in the list is technically one row of data with named columns,
        represented as a dictionary object
        For example:
        list = [
            {"col1":1234567, "col2":1234, "col3":123456, "col4":BLAH},
            {"col1":7654321, "col2":1234, "col3":123456, "col4":BLAH}
            ]
            """
        self.exec(sql,bindvars=bindvars)
#        if self.cursor.rowcount > 0: DOES NOT Work for Oracle until cursor.fetchone() or cursor.fetchall()
        columns = [i[0] for i in self.cursor.description]
        new_list = []
        for row in self.cursor:
            row_dict = dict()
            for col in columns:
                if upCase:
                    row_dict[col.upper()] = row[columns.index(col)]
                elif lowCase:
                    row_dict[col.lower()] = row[columns.index(col)]
                else:
                    row_dict[col] = row[columns.index(col)]

            new_list.append(row_dict)
        return new_list

    def nmget_dict_list(self, sql, bindvars=None, upCase=True, batchsize=100):
        mLst = []
        for bv in bindvars:
            mLst += self.nget_dict_list(self, sql, bindvars=bv, upCase=upCase, batchsize=batchsize)
        return(mLst)

    def mget_dict_list(self, sql, bindvars=None,upCase=True):
        new_list = []
        columns = []
        for bv in bindvars:
            self.exec(sql,bindvars=bv)
#            if self.cursor.rowcount > 0: DOES NOT Work for Oracle until cursor.fetchone() or cursor.fetchall()
            if len(columns)<1:
                columns = [i[0] for i in self.cursor.description]
            for row in self.cursor:
                row_dict = dict()
                for col in columns:
                    if upCase:
                        row_dict[col.upper()] = row[columns.index(col)]
                    else:
                        row_dict[col] = row[columns.index(col)]
                new_list.append(row_dict)
        return new_list

    def get_dataframe(self, sql, bindvars=None,upCase=False,lowCase=True):
        new_dataframe = pd.DataFrame.from_dict(self.get_dict_list(sql, bindvars=bindvars,upCase=upCase,lowCase=lowCase))
        new_dataframe.columns = [x.replace(' ','_').replace('\'','') for x in new_dataframe.columns]
        return new_dataframe

    def add_dataframe(self, df, sql, bindvars=None):
        for idx,row in df.iterrows():
            ndict = self.get_single_dict(sql,row)
        new_dataframe = pd.DataFCherame.from_dict(self.get_dict_list(sql, bindvars=None))
        new_dataframe.columns = [x.replace(' ','_').replace('\'','') for x in new_dataframe.columns]
        return ndict

    # -----------------------------------------------------------------------------------------------------
    # SQL Generation
    # -----------------------------------------------------------------------------------------------------
    def gen_sqltyping(self,var,varname=None):
        raise NotImplementedError("Each SqlConnector is responsible for its own gen_sqltyping method.")

    def fix_sqlvars(self,sql,bindvars):
        rx = re.compile(r":\w+")
        args = [m.group()[1:] for m in rx.finditer(sql)]
        newvars = {}
        for arg in args:
            newvars[arg] = bindvars[arg]
        return newvars

    def sqlquote(self,str):
        return("'"+str+"'")

    def sqlbindvar(self,str):
        return(":"+str)

    def gen_UpdateSQL(self,uTbl,uDict,uWhere,bindvars=False,includeNone=False):
        uSql = "Update " + uTbl + " Set "
        if uWhere:
            uWhere = " Where " + uWhere
        lSet = []
        for p in uDict:
            if uDict.get(p) or includeNone:
            #if uDict[p] is not None and excludeNone:
                if bindvars:
                    lSet.append(p + " = " + self.gen_sqltyping(uDict[p],p))
                else:
                    lSet.append(p + " = " + self.gen_sqltyping(uDict[p]))
        uSet = ', '.join(map(str, lSet))
        return(uSql + uSet + uWhere)

    # -----------------------------------------------------------------------------------------------------
    # Upload Procedures
    # -----------------------------------------------------------------------------------------------------
    def upload_entry(self,uplTable,uplWhereLst,uplDataDict,overwrite=False):
        retDict = {}

        wWhere = " and ".join([f"{w} = {self.sqlbindvar(w)}" for w in uplWhereLst])
        wDict = {k: v for k, v in uplDataDict.items() if k in uplWhereLst}

        cntSql = f"Select count(1) From {uplTable} Where {wWhere}"

        insCol = ", ".join([f"{w}" for w in uplWhereLst])
        insVal = ", ".join([f"{self.sqlbindvar(w)}" for w in uplWhereLst])
        insSQL = f"Insert into {uplTable} ({insCol}) Values ({insVal})"
        
        updSetLst = [k for k in uplDataDict if k not in uplWhereLst]
        updSet = ", ".join([f"{w} = {self.sqlbindvar(w)}" for w in updSetLst])
        updSQL = f"Update {uplTable} Set {updSet} Where {wWhere}"

        nCnt = self.nCount(cntSql, bindvars=wDict )
        retDict['Log'] = f"[uploadSQL] {uplTable} - "
        if nCnt == 0:
            logger.debug(f"[uploadSQL] {uplTable} - Insert {wDict} ")
            self.exec(insSQL , bindvars=wDict, commit=True)
            retDict['Log'] += "Insert "
            overwrite = True
        if overwrite:
            logger.debug(f"[uploadSQL] {uplTable} - Update {wDict} ")
            self.exec(updSQL , bindvars=uplDataDict, commit=True)
            retDict['Log'] += "Update "
        return(retDict)

    def upload_entry_list(self,uplTable,uplWhereLst,uplLst,overwrite=False):
        if len(uplLst) > 0:
            wWhere = " and ".join([f"{w} = {self.sqlbindvar(w)}" for w in uplWhereLst])
            cntSql = f"Select count(1) From {uplTable} Where {wWhere}"

            insCol = ", ".join([f"{w}" for w in uplWhereLst])
            insVal = ", ".join([f"{self.sqlbindvar(w)}" for w in uplWhereLst])
            insSQL = f"Insert into {uplTable} ({insCol}) Values ({insVal})"

            for uplDataDict in uplLst:
                wDict = {k: v for k, v in uplDataDict.items() if k in uplWhereLst}
                updSetLst = [k for k in uplDataDict if k not in uplWhereLst]

                updSet = ", ".join([f"{w} = {self.sqlbindvar(w)}" for w in updSetLst])
                updSQL = f"Update {uplTable} Set {updSet} Where {wWhere}"
                nCnt = self.nCount(cntSql, bindvars=wDict )
                if nCnt == 0:
                    self.exec(insSQL , bindvars=wDict, commit=True)
                    overwrite = True
                if overwrite:
                    self.exec(updSQL , bindvars=uplDataDict, commit=True)


