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
logName = "Upload_Collaborator"
#logFileName = os.path.join(djDir,"applog",f"x{logName}_{logTime:%Y%m%d_%H%M%S}.log")
logLevel = logging.INFO 

logger = logging.getLogger(__name__)
logging.basicConfig(
    format="[%(name)-20s] %(message)s ",
#    handlers=[logging.FileHandler(logFileName,mode='w'),logging.StreamHandler()],
    handlers=[logging.StreamHandler()],
    level=logLevel)
#-----------------------------------------------------------------------------


def get_oraCollaborator(test=0):
    from oraCastDB.oraCastDB import openCastDB

    usrSQL = "Select * From  Collaborator"
    grpSQL = "Select * From  Collaborator_Group"

    if test>0:
        usrSQL += f" Fetch First {test} Rows Only "
        grpSQL += f" Fetch First {test} Rows Only "

    CastDB = openCastDB()
    logger.info(f"[Collaborator] ... ")
    collabDict = {}
    #collabDict['User'] = pd.DataFrame(CastDB.get_dict_list(usrSQL))
    collabDict['Group'] = pd.DataFrame(CastDB.get_dict_list(grpSQL))

    uploadDir = 'C:/Code/zdjCode/adjCOADD/utilities/upload_data/Data'
    XlsFile = os.path.join(uploadDir,'CollaboratorData_v01.xlsx')
    if os.path.exists(XlsFile):

        # Organisation 
        XlsSheet = 'Organisation'
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        collabDict['Organisation'] = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        # User
        XlsSheet = 'User'
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        collabDict['User'] = pd.read_excel(XlsFile, sheet_name=XlsSheet)

        # User
        XlsSheet = 'Group'
        logger.info(f"[adjCOADD] Read {XlsFile}[{XlsSheet}] ")
        collabDict['Group'] = pd.read_excel(XlsFile, sheet_name=XlsSheet)

    logger.info(f"[Collaborators] {len(collabDict['User'])} ")
    logger.info(f"[Groups       ] {len(collabDict['Group'])} ")
    logger.info(f"[Organisation ] {len(collabDict['Organisation'])} ")
    CastDB.close()


    # logger.info(f"DF - Rename Columns {len(renameCol)}")
    # collabDict['User'].rename(columns=renameCol, inplace=True)
    # collabDict['Group'].rename(columns=renameCol, inplace=True)

    return(collabDict)

#-----------------------------------------------------------------------------
def main(prgArgs,djDir):

    sys.path.append(djDir)
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "adjcoadd.settings")
    django.setup()

    from apputil.models import Dictionary
    from applib.data.set_fielddata import set_model_arrayfields, set_dictFields, set_model_dicts
    from dcollab.models import Organisation, Collab_User, Collab_Group

    from django_countries import countries
    
    logger.info(f"Python         : {sys.version.split('|')[0]}")
    logger.info(f"Conda Env      : {os.environ['CONDA_DEFAULT_ENV']}")
    #logger.info(f"LogFile        : {logFileName}")

    logger.info(f"Django         : {django.__version__}")
    logger.info(f"Django Folder  : {djDir}")
    logger.info(f"Django Project : {os.environ['DJANGO_SETTINGS_MODULE']}")

   # Table -------------------------------------------------------------

    if prgArgs.table == "Collaborator" :

        Run = ['Group']

        print("--> oraCollaborator -----------------------------------------------------")
        CollabData = get_oraCollaborator(int(prgArgs.test))
        
        if "Organisation" in Run:
            print("-ORGANISATION --------------------------------------------------------------")
            print(CollabData["Organisation"].columns)
            OrganLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in CollabData["Organisation"].to_dict(orient='records')]
            nProcessed = 0
            for org in tqdm(OrganLst):
                nProcessed = nProcessed + 1
                #print(stock)

                newEntry = False
                djOrg = Organisation.get(None,org['organisation_name'],None)
                if djOrg is None:
                    newEntry = True
                    djOrg = Organisation()
                    djOrg.organisation_name = org['organisation_name']

                djOrg.organisation_code = org['organisation_code']
                djOrg.organisation_type = Dictionary.get(djOrg.DICTIONARY_FIELDS["organisation_type"],org['organisation_type'])

                for code, name in list(countries):
                    if name == org['country']:
                        djOrg.country = code
                if djOrg.country is None:
                    print(f"No Country Code for {org['country']}")

                if prgArgs.upload:
                    if prgArgs.overwrite or newEntry:
                        print(f"{repr(djOrg)} {newEntry}")
                        djOrg.save()

        if "User" in Run:
            print("-USER -----------------------------------------------------------------")
            print(CollabData["User"].columns)
            UserLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in CollabData["User"].to_dict(orient='records')]
            cpyFields = ['ora_user_id', 'ora_group_id', 'title', 'first_name', 'last_name',
                        'position', 'department', 'postal_address',
                        'country', 'phone',  'email1', 'email2', 'active_email']
            
            nProcessed = 0
            for user in tqdm(UserLst):
                nProcessed = nProcessed + 1
                #print(stock)

                newEntry = False
                if 'email1' in user:
                    djUser = Collab_User.get(None,user['email1'],None,None)
                else:
                    djUser = Collab_User.get(None,None,user['first_name'],user['last_name'])
                if djUser is None:
                    newEntry = True
                    djUser = Collab_User()

                for f in cpyFields:
                    if f in user:
                        setattr(djUser, f, user[f])

                djOrg = Organisation.get(None,user['organisation'])
                if djOrg is None:
                    print(f"No Organisation for {user['organisation']}")
                else:
                    djUser.organisation_id = djOrg

                for code, name in list(countries):
                    if name == user['country'].strip():
                        djUser.country = code
                if djUser.country is None:
                    print(f"No Country Code for {user['country']}")

                if prgArgs.upload:
                    if prgArgs.overwrite or newEntry:
                        print(f"{repr(djUser)} {newEntry} ({djUser.country})")
                        djUser.save()



        if "Group" in Run:
            print("-GROUP -----------------------------------------------------------------")
            print(CollabData["Group"].columns)
            GroupLst = [{k:v for k,v in m.items() if pd.notnull(v)} for m in CollabData["Group"].to_dict(orient='records')]

            cpyFields = ['ora_group_id', 'ora_pi_id', 'group_code',
                        'department', 'postal_address',
                        'country', 'email', 'mta_document']


            nProcessed = 0
            for group in tqdm(GroupLst):
                nProcessed = nProcessed + 1

                newEntry = False
                if 'group_code' in group:
                    djGroup = Collab_Group.get(None,group['group_code'],None,None)
                else:
                    djGroup = Collab_Group.get(None,None,group['pi_user_id'],group['organisation_id'])
                if djGroup is None:
                    newEntry = True
                    djGroup = Collab_Group()

                for f in cpyFields:
                    if f in group:
                        setattr(djGroup, f, group[f])

                djGroup.mta_status = Dictionary.get(djGroup.DICTIONARY_FIELDS["mta_status"],group['mta_status'])

                djOrg = Organisation.get(None,group['organisation'])
                if djOrg is None:
                    print(f"No Organisation for {group['organisation']}")
                else:
                    djGroup.organisation_id = djOrg

                for code, name in list(countries):
                    if name == group['country'].strip():
                        djGroup.country = code
                if djGroup.country is None:
                    print(f"No Country Code for {group['country']}")

                if prgArgs.upload:
                    if prgArgs.overwrite or newEntry:
                        #print(f"{repr(djGroup)} {newEntry}")
                        djGroup.save()


#==============================================================================
if __name__ == "__main__":

    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")


    # ArgParser -------------------------------------------------------------
    prgParser = configargparse.ArgumentParser(prog='upload_Django_Data', 
                                description="Uploading data to adjCOADD from Oracle/Excel/CSV")
    prgParser.add_argument("-t",default=None,required=True, dest="table", action='store', help="Table to upload [CompoundID]")
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

    prgArgs = prgParser.parse_args()

    # Django -------------------------------------------------------------
    if prgArgs.django == 'Meran':
        djDir = "D:/Code/zdjCode/adjCOADD"
    #   uploadDir = "C:/Code/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    #   orgdbDir = "C:/Users/uqjzuegg/The University of Queensland/IMB CO-ADD - OrgDB"
    elif prgArgs.django == 'Work':
        djDir = "/home/uqjzuegg/xhome/Code/zdjCode/adjCOADD"
    #     uploadDir = "C:/Data/A02_WorkDB/03_Django/adjCOADD/utilities/upload_data/Data"
    elif prgArgs.django == 'Laptop':
        djDir = "C:/Code/zdjCode/adjCOADD"
    #     uploadDir = "/home/uqjzuegg/DeepMicroB/Code/Python/Django/adjCOADD/utilities/upload_data/Data"
    else:
        djDir = None

    if djDir:
        main(prgArgs,djDir)
        print("-------------------------------------------------------------------")

#==============================================================================
