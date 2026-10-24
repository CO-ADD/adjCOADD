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
def set_model_fkeys(djObj, rowDict, dict_FKeys, **kwargs):
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
                    setattr(djObj,f,_obj)
                else:
                    if verbose >0:
                        logger.warning(f" [set_fkeys] {f} = '{rowDict[f]}' not found in [{dict_FKeys[f].__name__}]" )
                    if valLog:
                        valLog.add_error(f"FKey {f} not found",rowDict[f],f"Not found in [{dict_FKeys[f].__name__}]")
                        
                    setattr(djObj,f,None)
                    valid = False
    return valid           
        
#------------------------------------------------------------------------------------
def set_model_fields(djObj,rowDict,list_Fields,**kwargs):
    #
    # list_Fields = ['fieldname1','fielname2',...,'filenameN'] 
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    for e in list_Fields:
        if e in rowDict:
            if pd.notnull(rowDict[e]):
                setattr(djObj,e,rowDict[e])

#------------------------------------------------------------------------------------
def set_model_dicts(djObj,rowDict,list_Dicts,**kwargs):
    #
    # list_Dicts = ['fieldname1','fielname2',...,'filenameN'] 
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
        
    for d in list_Dicts:
        #print(f" [set_model_dicts] - {d}")
        if d in rowDict:
            #print(f" [set_model_dicts] - {d} {rowDict[d]} -> {djObj.DICTIONARY_FIELDS}")
            if pd.notnull(rowDict[d]):
                if d in djObj.DICTIONARY_FIELDS:
                    setattr(djObj,d,Dictionary.get(djObj.DICTIONARY_FIELDS[d],rowDict[d]))
    
# --------------------------------------------------------------------------------
    
    
    
#------------------------------------------------------------------------------------
def set_model_arrayfields(djObj,rowDict, dict_Arrays,**kwargs):
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
                setattr(djObj,f, _array)
                
        elif isinstance(dict_Arrays[f],list):
            for l in dict_Arrays[f]:
                if l in rowDict:
                    if pd.notnull(rowDict[l]):
                        _array.append(rowDict[l])
        if _array:
            setattr(djObj,f,_array)

#------------------------------------------------------------------------------------
def set_model_dictarrayfields(djObj,rowDict,dict_ArrayDicts,**kwargs):
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
                if f in djObj.DICTIONARY_FIELDS:
                    _dict_list = strList_to_List(rowDict[dict_ArrayDicts[f]])
                
        elif isinstance(dict_ArrayDicts[f],list):
            for l in dict_ArrayDicts[f]:
                if l in rowDict:
                    if pd.notnull(rowDict[l]):
                        _dict_list.append(rowDict[l])
        
        for l in _dict_list:
            _d = Dictionary.get(djObj.DICTIONARY_FIELDS[f],l)
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
                    valLog.add_error(f'Wrong {djObj.DICTIONARY_FIELDS[f]}',f"{l}",f"[{f}]",f"Correct the value {_help}")
            
        if len(_ret_list)>0:
            setattr(djObj,f,_ret_list)

#------------------------------------------------------------------------------------
def set_model_fkeyarrayfields(djObj,rowDict,dict_ArrayFKeys,**kwargs):
    #
    #  dict_ArrayFKeys = {'fieldname_lst' : {'model': Model, 'fields': ['input1','input2',...,'inputN'] } }
    #   inputN - checked if in Dictionary
    #
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)

    for f in dict_ArrayFKeys:
        #print("arrFields",f)
        _array = []
        _fields = dict_ArrayFKeys[f]['fields']
        _model = dict_ArrayFKeys[f]['model']
        
        if isinstance(_fields,str):
            if pd.notnull(rowDict[_fields]):
                _array = strList_to_List(str(rowDict[_fields]))
                setattr(djObj,f,_array)
                
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
                        if verbose >0:
                            logger.warning(f" [set_fkeyarrays] {f} = '{_a}' not found in [{_model.__name__}]" )
                        if valLog:
                            valLog.add_error(f"Missing {f}",_a,f"Not found in [{_model.__name__}]","Correct MotherPlate_ID or Upload MotherPlate")
                        _valid = False
            if _valid:
                setattr(djObj,f,_array)


#------------------------------------------------------------------------------------
def set_model_from_dict(djObj,row,
                        list_Fields=[], 
                        dict_Arrays={}, 
                        list_Dicts=[],
                        dict_FKeys={},
                        dict_FKeyArrays = {}, 
                        **kwargs):
    
    valLog = kwargs.get('valLog',None)
    verbose = kwargs.get('verbose',0)
    validStatus = True

    if len(list_Fields)>0:
        set_model_fields(djObj,row,list_Fields)
    if len(dict_Arrays)>0:
        set_model_arrayfields(djObj,row,dict_Arrays)     
    if len(list_Dicts)>0:
        set_model_dicts(djObj,row,list_Dicts,valLog=valLog)
    if len(dict_FKeys)>0:
        set_model_fkeys(djObj,row,dict_FKeys,valLog=valLog)
    if len(dict_FKeyArrays)>0:
        set_model_fkeyarrayfields(djObj,row,dict_FKeyArrays,valLog=valLog)
        
    djObj.set_defaults_model()
    validDict = djObj.validate_fields()
    if validDict:
        validStatus = False
        for k in validDict:
            if valLog:    
                valLog.add_log('Warning','',k,validDict[k],'-')
            if verbose >0: 
                logger.warning(f"{k} - {validDict[k]}")
    return(validStatus)

#------------------------------------------------------------------------------------
def rename_fkeyarray(djObj, oldPK, newPK, fkeyArray):
    #
    # fklModel - Model with ForeignKeyList
    # oldPK - old Primary Key
    # newPK - new Primary Key
    #
    #print(f" [rename_fkeyarray] {djObj.__class__.__name__}.{fkeyArray} {oldPK} -> {newPK} ")
    
    _fka = getattr(djObj,fkeyArray)
    if oldPK in _fka:
        idx = _fka.index(oldPK)
        _fka[idx] = newPK
        setattr(djObj,fkeyArray,_fka)