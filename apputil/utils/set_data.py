#
import pandas as pd
from apputil.utils.data import strList_to_List
from apputil.models import Dictionary
#

# dictList = ['project_name','project_comment',
#             ...
#             ]
# arrDict = {'screen_status': 'screen_status',                  # strList->List
#             'ora_contact_ids':['CONTACT_A_ID','CONTACT_B_ID'] # append.List
#            }
# dictFields = ['project_type','provided_container','stock_conc_unit',] # using Choice_Dictionary


def set_arrayFields(djModel,rowDict, arrDict):
    for f in arrDict:
        #print("arrFields",f)
        if isinstance(arrDict[f],str):
            if pd.notnull(rowDict[arrDict[f]]):
                setattr(djModel,f,strList_to_List(rowDict[arrDict[f]]))
                
        elif isinstance(arrDict[f],list):
            _list = []
            for l in arrDict[f]:
                if pd.notnull(rowDict[l]):
                    _list.append(rowDict[l])
            setattr(djModel,f,_list)
        
def set_dictFields(djModel,rowDict,dictList):
    for e in dictList:
        #print("fromDict",e)
        if pd.notnull(rowDict[e]):
            setattr(djModel,e,rowDict[e])

def set_Dictionaries(djModel,rowDict,dictFields):
    for d in dictFields:
        #print("fromDict",d)
        if pd.notnull(rowDict[d]):
            if d in djModel.Choice_Dictionary:
                setattr(djModel,d,Dictionary.get(djModel.Choice_Dictionary[d],rowDict[d]))            


# if upload:
#    if overwrite:
#       save
#    elif newentry:
#       save

def set_fields_fromrow(djModel,row,dictList=[], arrDict=[], dictFields=[],valLog=None):
    if len(dictList)>0:
        set_dictFields(djModel,row,dictList)
    if len(arrDict)>0:
        set_arrayFields(djModel,row,arrDict)     
    if len(dictFields)>0:
        set_Dictionaries(djModel,row,dictFields)
        
    djModel.clean_Fields()
    validDict = djModel.validate()
    if validDict:
        validStatus = False
        for k in validDict:
            if valLog:    
                valLog.add_log('Warning','',k,validDict[k],'-')
            else: 
                print('Warning',k,validDict[k],'-')