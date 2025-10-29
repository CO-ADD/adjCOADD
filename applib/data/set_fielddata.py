#
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

from applib.data.str_lists import strList_to_List
from apputil.models import Dictionary
#

# dictList = ['project_name','project_comment',
#             ...
#             ]
# arrDict = {'screen_status': 'screen_status',                  # strList->List
#             'ora_contact_ids':['CONTACT_A_ID','CONTACT_B_ID'] # append.List
#            }
# dictFields = ['project_type','provided_container','stock_conc_unit',] # using DICTIONARY_FIELDS


#------------------------------------------------------------------------------------
def set_model_fkeys(djModel, rowDict, dict_FKeys, **kwargs):
    #
    # dict_FKeys = {'field_name': FKey Model}
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    valid = True
    for f in dict_FKeys:
        if f in rowDict:
            if rowDict[f] != '-' and rowDict[f] is not None:
                #print(f"{rowDict[f]} {type(rowDict[f])}")
                _obj = dict_FKeys[f].get(rowDict[f])
                if _obj is not None:
                    setattr(djModel,f,_obj)
                else:
                    logger.warning(f" [set_fkeys] {f} = '{rowDict[f]}' not found in [{dict_FKeys[f].__name__}]" )
                    if valLog:
                        valLog.add_error(f"FKey {f} not found",rowDict[f],"Not found in [{dict_FKeys[f].__name__}]")
                        
                    setattr(djModel,f,None)
                    valid = False
    return valid           
        
#------------------------------------------------------------------------------------
def set_model_fields(djModel,rowDict,list_Fields,**kwargs):
    #
    # list_Fields = ['fieldname1','fielname2',...,'filenameN'] 
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    for e in list_Fields:
        if e in rowDict:
            if pd.notnull(rowDict[e]):
                setattr(djModel,e,rowDict[e])

#------------------------------------------------------------------------------------
def set_model_dicts(djModel,rowDict,list_Dicts,**kwargs):
    #
    # list_Dicts = ['fieldname1','fielname2',...,'filenameN'] 
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
        
    for d in list_Dicts:
        #print(f" [set_model_dicts] - {d}")
        if d in rowDict:
            #print(f" [set_model_dicts] - {d} {rowDict[d]} -> {djModel.DICTIONARY_FIELDS}")
            if pd.notnull(rowDict[d]):
                if d in djModel.DICTIONARY_FIELDS:
                    setattr(djModel,d,Dictionary.get(djModel.DICTIONARY_FIELDS[d],rowDict[d]))
    
# --------------------------------------------------------------------------------
    
    
    
#------------------------------------------------------------------------------------
def set_model_arrayfields(djModel,rowDict, dict_Arrays,**kwargs):
    #
    # dict_Arrays = {'fieldname_lst' : ['input1','input2',...,'inputn'] }
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    for f in dict_Arrays:
        _array = []
        if isinstance(dict_Arrays[f],str):
            if pd.notnull(rowDict[dict_Arrays[f]]):
                _array = strList_to_List(str(rowDict[dict_Arrays[f]]))
                setattr(djModel,f, _array)
                
        elif isinstance(dict_Arrays[f],list):
            for l in dict_Arrays[f]:
                if l in rowDict:
                    if pd.notnull(rowDict[l]):
                        _array.append(rowDict[l])
        if _array:
            setattr(djModel,f,_array)

#------------------------------------------------------------------------------------
def set_model_dictarrayfields(djModel,rowDict,dict_ArrayDicts,**kwargs):
    #
    # arrDict = {'fieldname_lst' : ['input1','input2',...,'inputN'] }
    #   inputN - checked if in Dictionary
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    for f in dict_ArrayDicts:
        _dict_list = []
        _ret_list  = []
        if isinstance(dict_ArrayDicts[f],str):
            if pd.notnull(rowDict[dict_ArrayDicts[f]]):
                if f in djModel.DICTIONARY_FIELDS:
                    _dict_list = strList_to_List(rowDict[dict_ArrayDicts[f]])
                
        elif isinstance(dict_ArrayDicts[f],list):
            for l in dict_ArrayDicts[f]:
                if l in rowDict:
                    if pd.notnull(rowDict[l]):
                        _dict_list.append(rowDict[l])
        
        for l in _dict_list:
            _d = Dictionary.get(djModel.DICTIONARY_FIELDS[f],l)
            if _d:
                _ret_list.append(str(_d))
            else:
                if valLog:
                    if 'conc' in f:
                        _help = " [µM -> uM]"
                    elif 'amount' in f:
                        _help = " [µg -> ug]"
                    elif 'volume' in f:
                        _help = " [µl -> uL]"
                    else:
                        _help = ""
                    valLog.add_error(f'Wrong {djModel.DICTIONARY_FIELDS[f]}',f"{l}",f"[{f}]",f"Correct the value {_help}")
            
        if len(_ret_list)>0:
            setattr(djModel,f,_ret_list)

#------------------------------------------------------------------------------------
def set_model_fkeyarrayfields(djModel,rowDict,dict_ArrayFKeys):
    #
    #  dict_ArrayFKeys = {'fieldname_lst' : {'model': Model, 'fields': ['input1','input2',...,'inputN'] } }
    #   inputN - checked if in Dictionary
    #

    for f in dict_ArrayFKeys:
        #print("arrFields",f)
        _array = []
        _fields = dict_ArrayFKeys[f]['fields']
        _model = dict_ArrayFKeys[f]['model']
        
        if isinstance(_fields,str):
            if pd.notnull(rowDict[_fields]):
                _array = strList_to_List(str(rowDict[_fields]))
                setattr(djModel,f,_array)
                
        elif isinstance(_fields,list):
            for l in _fields:
                if l in rowDict:
                    if pd.notnull(rowDict[l]):
                        _array.append(rowDict[l])
        if _array:
            _valid = True
            for _a in _array:
                if _a != '-':
                    _obj = _model.get(_a)
                    if _obj is None:
                        logger.warning(f" [set_fkeyarrays] {f} = '{_a}' not found in [{_model.__name__}]" )
                        _valid = False
            if _valid:
                setattr(djModel,f,_array)


#------------------------------------------------------------------------------------
def set_model_from_dict(djModel,row,
                        list_Fields=[], 
                        dict_Arrays={}, 
                        list_Dicts=[],
                        dict_FKeys={},
                        dict_FKeyArrays = {}, 
                        valLog=None):
    validStatus = True

    if len(list_Fields)>0:
        set_model_fields(djModel,row,list_Fields)
    if len(dict_Arrays)>0:
        set_model_arrayfields(djModel,row,dict_Arrays)     
    if len(list_Dicts)>0:
        set_model_dicts(djModel,row,list_Dicts)
    if len(dict_FKeys)>0:
        set_model_fkeys(djModel,row,dict_FKeys)
    if len(dict_FKeyArrays)>0:
        set_model_fkeyarrayfields(djModel,row,dict_FKeyArrays)
        
    djModel.set_defaults_model()
    validDict = djModel.validate_fields()
    if validDict:
        validStatus = False
        for k in validDict:
            if valLog:    
                valLog.add_log('Warning','',k,validDict[k],'-')
            else: 
                logger.warning(f"{k} - {validDict[k]}")
    return(validStatus)
