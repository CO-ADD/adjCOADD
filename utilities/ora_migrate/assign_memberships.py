#
#
#
import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
import configargparse
from pathlib import Path

from tqdm import tqdm
# from zUtils import zData

import django
#from djCOADD import djOrgDB
from oraCastDB import oraCastDB
#-----------------------------------------------------------------------------

# Logger ----------------------------------------------------------------
import logging
logTime= datetime.datetime.now()
logName = "Assign_Membership"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    django.setup()

    from apputil.models import Dictionary
    #from applib.data.set_fielddata import set_model_arrayfields, set_model_dicts
    from dcollab.models import Organisation, Collab_User, Collab_Group, Collab_Membership
    from dsample.models import Project, Project_Membership

    from django_countries import countries
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "Group" :

        qryGroup = Collab_Group.objects.all()
        for grp in qryGroup:
            if grp.ora_pi_id:
                pi = Collab_User.objects.filter(ora_user_id=grp.ora_pi_id)
                if pi:
                    _m, _created = Collab_Membership.objects.get_or_create(user_id=pi[0],group_id=grp)
                    _m.role='LI'
                    _m.save()
                    print(repr(_m))

    elif prgArgs.table == "User" :
        qryUser = Collab_User.objects.all()
        for usr in qryUser:
            if usr.ora_group_id:
                grp = Collab_Group.objects.filter(ora_group_id=usr.ora_group_id)
                if grp:
                    _m, _created = Collab_Membership.objects.get_or_create(user_id=usr,group_id=grp[0])
                    if _m != 'LI':
                        _m.role = 'M'
                        _m.save
                    print(repr(_m))

    elif prgArgs.table == "Project" :
        qryProject = Project.objects.all()
        for prg in qryProject:
            n = 1
            if prg.ora_contact_ids:
                for contact in prg.ora_contact_ids:
                    _c = Collab_User.objects.filter(ora_user_id=contact)
                    if _c:
                        print(f" {prg} {contact} {_c[0]} {n}")
                        _m, _created = Project_Membership.objects.get_or_create(user_id=_c[0],project_id=prg)
                        if n == 1:
                            _m.role = 'PC'
                            _m.save()
                        n +=1 

#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [Group/Project]")
    prgParser.add_argument("--upload",default=False,required=False, dest="upload", action='store_true', help="Upload data to dj Database")
    prgParser.add_argument("--overwrite",default=False,required=False, dest="overwrite", action='store_true', help="Overwrite existing data")
    prgParser.add_argument("--user",default='J.Zuegg',required=False, dest="appuser", action='store', help="AppUser to Upload data")
    prgParser.add_argument("--test",default=0,required=False, dest="test", action='store', help="Number of entries to test")

#    prgParser.add_argument("-d","--directory",default=None,required=False, dest="directory", action='store', help="Directory or Folder to parse")
#    prgParser.add_argument("-f","--file",default=None,required=False, dest="file", action='store', help="Single File to parse")
#    prgParser.add_argument("--db",default='Local',required=False, dest="database", action='store', help="Database [Local/Work/WorkLinux]")
#    prgParser.add_argument("-r","--runid",default=None,required=False, dest="runid", action='store', help="Antibiogram RunID")

    prgParser.add_argument("--django",default='Local',required=False, dest="django", action='store', help="Django configuration [Meran/Laptop/Work]")
    prgParser.add_argument("-c","--config",type=Path,is_config_file=True,help="Path to a configuration file ",)

    try:
        prgArgs = prgParser.parse_args()
    except:
        prgParser.print_help()
        sys.exit(0)

    from zDjango.djUtils import init_django_dir

    # Django -------------------------------------------------------------
    djDir = init_django_dir(prgArgs,"adjCOADD")
    if djDir:
        print(djDir)
        main(prgArgs,djDir)

#==============================================================================
