import os,sys,io
import zipfile
#from datetime import time, date
import csv
import numpy as np
import pandas as pd

from apputil.utils.data import listFolders
#-----------------------------------------------------------------------------
def parse_FastQC(AssemblyFolder,OrgBID,RunID):
#-----------------------------------------------------------------------------
    BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
    FastqcDir = os.path.join(BatchDir,"trim","paired","fastqc")
    outLst = []

    P1Zip = f"{OrgBID}_{RunID}_P1.zip"
    with zipfile.ZipFile(os.path.join(FastqcDir,P1Zip)) as z:
        for filename in z.namelist():
            if "summary.txt" in filename:
                pDict ={}
                pDict['Seq'] = 'P1'
                with z.open(filename,'r') as f:
                    reader = csv.reader(io.TextIOWrapper(f),delimiter='\t')
                    for row in reader:
                        pDict[row[1]] = row[0]
                    outLst.append(pDict)

    P2Zip = f"{OrgBID}_{RunID}_P2.zip"
    with zipfile.ZipFile(os.path.join(FastqcDir,P2Zip)) as z:
        for filename in z.namelist():
            if "summary.txt" in filename:
                pDict ={}
                pDict['Seq'] = 'P2'
                with z.open(filename,'r') as f:
                    reader = csv.reader(io.TextIOWrapper(f),delimiter='\t')
                    for row in reader:
                        pDict[row[1]] = row[0]
                    outLst.append(pDict)

    return(outLst)
#-----------------------------------------------------------------------------
def parse_CheckM(AssemblyFolder,OrgBID,RunID,outType="scaffolds"):
#-----------------------------------------------------------------------------
    BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
    CheckmDir = os.path.join(BatchDir,"assembly","spades","checkM","lineage","storage")
    outDict = {}
    if os.path.exists(CheckmDir):
        with open(os.path.join(CheckmDir,"bin_stats_ext.tsv")) as file:
            tsv_file = csv.reader(file,delimiter="\t")
            binDict = {}
            for line in tsv_file:
                binDict[line[0]]=eval(line[1])
        
            outDict["marker_lineage"] = binDict[outType]["marker lineage"]
            outDict["n_genomes"] = binDict[outType]["# genomes"]
            outDict["n_predit_genes"] = binDict[outType]["# predicted genes"]
            outDict["n_markers"] = binDict[outType]["# markers"]
            outDict["n_marker_sets"] = binDict[outType]["# marker sets"]
            outDict["n_scaffolds"] = binDict[outType]["# scaffolds"]
            outDict["genome_size"] = binDict[outType]["Genome size"]
            outDict["coding_density"] = binDict[outType]["Coding density"]
            outDict["completeness"] = binDict[outType]["Completeness"]
            outDict["contamination"] = binDict[outType]["Contamination"]
            outDict["genome_size"] = binDict[outType]["Genome size"]
            outDict["gc"] = binDict[outType]["GC"]
            outDict["gc_std"] = binDict[outType]["GC std"]
            outDict["n_ambig_bases"] = binDict[outType]["# ambiguous bases"]
            outDict["longest_scaffold"] = binDict[outType]["Longest scaffold"]
            outDict["mean_scaffols"] = binDict[outType]["Mean scaffold length"]
            outDict["N50"] = binDict[outType]["N50 (scaffolds)"]
            outDict["trans_table"] = binDict[outType]["Translation table"]
    return(outDict)


#-----------------------------------------------------------------------------
def parse_Kraken(AssemblyFolder,OrgBID,RunID,outType="S",inType="fasta",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
    KrakenDir = os.path.join(BatchDir,"kraken")
    KrakenF = f"{OrgBID}_{RunID}_{inType}.report"
    outLst = []
    if os.path.exists(KrakenDir):
        with open(os.path.join(KrakenDir,KrakenF)) as file:
            tsv_file = csv.reader(file,delimiter="\t")
            for line in tsv_file:
                if line[3] == outType:
                    pctSeq = float(line[0])
                    if pctSeq >= pctCutOff:
                        org_name = line[5].strip()
                        outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
    return(outLst)

#-----------------------------------------------------------------------------
def parse_Abricate(AssemblyFolder,OrgBID,RunID,inType="fasta",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
    AbricateDir = os.path.join(BatchDir,"abricate")
    AbricateF = f"{OrgBID}_{RunID}_{inType}.out"
    outLst = []
    if os.path.exists(AbricateDir):
        with open(os.path.join(AbricateDir,AbricateF)) as file:
            tsv_file = csv.reader(file,delimiter="\t")
            for line in tsv_file:
                if line[3] == outType:
                    pctSeq = float(line[0])
                    if pctSeq >= pctCutOff:
                        org_name = line[5].strip()
                        outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
    return(outLst)


#-----------------------------------------------------------------------------
def parse_MLST(AssemblyFolder,OrgBID,RunID,outType="S",inType="fastq",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
    KrakenDir = os.path.join(BatchDir,"mlst")
    KrakenF = f"{OrgBID}_{RunID}_{inType}.report"
    outLst = []
    if os.path.exists(KrakenDir):
        with open(os.path.join(KrakenDir,KrakenF)) as file:
            tsv_file = csv.reader(file,delimiter="\t")
            for line in tsv_file:
                if line[3] == outType:
                    pctSeq = float(line[0])
                    if pctSeq >= pctCutOff:
                        org_name = line[5].strip()
                        outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
    return(outLst)




#-----------------------------------------------------------------------------
def parse_WGS_COADD(zAssemblyBase):
#-----------------------------------------------------------------------------

    lstFastQC = []
    lstCheckM = []
    lstKraken = []

    for dirHead in listFolders(zAssemblyBase):
        zAssemblyFolder = os.path.join(zAssemblyBase,dirHead)
        for BatchID in listFolders(zAssemblyFolder):            
            for RunID in listFolders(os.path.join(zAssemblyFolder,BatchID)):

                dirAss = os.path.join(zAssemblyFolder,BatchID,RunID)
                print(f"[WGS-COADD] {dirHead} {BatchID} {RunID} --> {dirAss}")

                # FastQC
                lFastQC=parse_FastQC(zAssemblyFolder,BatchID,RunID)
                sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID}
                for row in lFastQC:
                    lstFastQC.append(dict(sDict,**row))

                # CheckM
                dCheckM = parse_CheckM(zAssemblyFolder,BatchID,RunID)
                sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID}
                lstCheckM.append(dict(sDict,**dCheckM))

                # Kraken
                lKraken = parse_Kraken(zAssemblyFolder,BatchID,RunID,outType="S",inType="fasta")
                sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}
                
                # Abricate
                lAbricate = parse_abricate(zAssemblyFolder,BatchID,RunID,outType="S",inType="fasta")
                sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}

                # MLST
                lMlst = parse_Kraken(zAssemblyFolder,BatchID,RunID,outType="S",inType="fastq")
                sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}

                idLst = []
                for i in range(len(lKraken)):
                    v = lKraken[i]
                    idLst.append(f"{v['org_name']} ({v['pct']:.1f} pct)")
                sDict['ID_Organisms'] = idLst
                lstKraken.append(sDict)

        #         break
        #     break
        # break

    return(lstFastQC,lstCheckM,lstKraken)