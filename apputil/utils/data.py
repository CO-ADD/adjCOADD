#
# General utilities for processing List, Dictionaries
#


#-----------------------------------------------------------------------------------
def split_StrList(strList,sep=";"):
#-----------------------------------------------------------------------------------
    if strList:
        retLst = str(strList).split(sep)
        for i in range(len(retLst)):
            retLst[i] = retLst[i].strip()
    else:
        retLst = None
    return(retLst)

