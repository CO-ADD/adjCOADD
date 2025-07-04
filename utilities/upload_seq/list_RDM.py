import os, sys
import datetime
import csv
import pandas as pd
import numpy as np
from pathlib import Path

# ----------------------------------------------------
RDM_DIR = 'M:/Sequence/WGS'
RDM_DIR = '/home/uqjzuegg/RDM/MICROORGDB/Sequence/WGS'

AssemblyBase = os.path.join(RDM_DIR,'03_Fasta')
CSV_FiLE = 'RDM_FastA_List.csv'
# ----------------------------------------------------
def listFolders(Path):
    return([ name for name in os.listdir(Path) if os.path.isdir(os.path.join(Path, name)) ])

def listFiles(Path,Extension=None):
    files = [ name for name in os.listdir(Path) if os.path.isfile(os.path.join(Path, name)) ]
    if Extension:
        return([ file for file in files if file.endswith(Extension)])
    else:
        return(files)


def split_BatchID_RunID(batch_run_id):
    arrStr = batch_run_id.split("_")
    batchID = '_'.join(arrStr[0:3])
    runID = '_'.join(arrStr[3:])
    return batchID, runID


#===========================================================================
lst_WGS = []

for subDir in listFolders(AssemblyBase):
    for SeqID in listFolders(os.path.join(AssemblyBase,subDir)):
        # Ignore Folder with _ at the beginning
        if SeqID[0] != '_':
            seqDir = os.path.join(AssemblyBase,subDir,SeqID)
            OrgBatchID, RunID = split_BatchID_RunID(SeqID)
            nFasta = len(listFiles(seqDir,Extension=".fasta"))
            _WGS = {'ORGBATCH_ID':OrgBatchID, 'SEQRUN_ID':RunID,'N_FASTA':nFasta,'SUB_DIR':subDir,'SEQ_DIR':seqDir}
            lst_WGS.append(_WGS)

df_WGS = pd.DataFrame(lst_WGS)
df_WGS.to_csv(CSV_FiLE)
#===========================================================================
