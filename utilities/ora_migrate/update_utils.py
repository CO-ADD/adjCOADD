#
#
#
from dsample.models import Convert_ProjectID, Convert_CompoundID

#-------------------------------------------------------------------
def convert_castdb_compoundid_from_dj(djID):
    validStatus = True
    if djID is not None:
        if 'MCC' in djID:
            djID = djID.replace('MCC','MCC_').replace("_",":")
        else:
            try:
                oraID = Convert_CompoundID.objects.get(compound_id = oraID).ora_compound_id
            except:
                oraID = None
                validStatus = False
    else:
        oraID = None
        validStatus = False
    
    return(oraID,validStatus)
#-------------------------------------------------------------------
def convert_castdb_compoundid_from_ora(oraID):
    validStatus = True
    if oraID is not None:
        if 'MCC_' in oraID:
            djID = oraID.replace('MCC_','MCC').replace(":","_")
        else:
            try:
                djID = Convert_CompoundID.objects.get(ora_compound_id = oraID).compound_id
            except:
                djID = None
                validStatus = False
    else:
        djID = None
        validStatus = False
    
    return(djID,validStatus)