
import os
import pandas as pd
#from django.apps import apps
from django.core.validators import validate_email
from django.core.exceptions import ValidationError
from dcollab.models import Collab_Group, Collab_User, Organisation
from dsample.models import Project_Membership

# #-----------------------------------------------------------------------------
# def Upload_Project_Collab(CollabDict,  upload=False, overwrite=False, **kwargs):
# # --------------------------------------------------------------------------------

#     valLog = kwargs.get('valLog',None)
#     verbose = kwargs.get('verbose',0)
    
#     for key in CollabDict:
#         djUsr = Collab_User.get(None,CollabDict[key]['email'])
#         djOrg = Organisation.get_bysimilarity(CollabDict[key]['organisation'])
        
#         if djOrg is None:
#             djOrg = Organisation()
#             djOrg.organisation_name = CollabDict[key]['organisation']
#             valLog.add_warning("New Organisation",CollabDict[key]['organisation'])
#         else:
#             valLog.add_info("Existing Organisation",f"{djOrg}")
        
#         if djUsr is None:
#             djUsr = Collab_User()   
#             djUsr.email = CollabDict[key]['email'] 
#             djUsr.first_name = CollabDict[key]['first_name'] 
#             djUsr.last_name = CollabDict[key]['last_name']
#             djUsr.organisation_id = djOrg 
#             valLog.add_warning("New Collaborator",f"{djUsr}")
#         else:
#             valLog.add_info("Existing Collaborator",f"{djUsr}")
            
#         if 'LI' == key:
#             djGrp = Collab_Group.get(None,Code=None, PI_ID=djUsr.user_id)
#             if djGrp is None:
#                 djGrp = Collab_Group()
#                 djGrp.group_code = f"{CollabDict[key]['last_name']}{CollabDict[key]['first_name'][0]}_{djOrg.organisation_code}"
#                 djGrp.organisation_id = djOrg
#                 valLog.add_warning("New Group",f"{djGrp}")
#             else:
#                 valLog.add_info("Existing Group",f"{djGrp}")
    