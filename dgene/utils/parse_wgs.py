import os,sys,io
import zipfile
#from datetime import time, date
import csv
import numpy as np
import pandas as pd

from apputil.utils.data import listFolders
#-----------------------------------------------------------------------------
def split_BatchID_RunID(batch_run_id):
#-----------------------------------------------------------------------------
    arrStr = batch_run_id.split("_")
    batchID = '_'.join(arrStr[0:3])
    runID = '_'.join(arrStr[3:])
    return batchID, runID

#-----------------------------------------------------------------------------
def get_FastQC_Info(AssemblyFolder,OrgBID,RunID):
#-----------------------------------------------------------------------------

    # print(f"[fastQC] {OrgBID} {RunID} ")

    retLst = []

    #-- FastQC ------------------------------------------------------------------------
    FastqcInfo = os.path.join(AssemblyFolder,"trim","fastqc","info_fastqc.csv ")
    if os.path.isfile(FastqcInfo):
        with open(FastqcInfo, 'r') as file:
            csv2dict = csv.DictReader(file)
            for fdict in list(csv2dict):
                fastqcDict = fdict
                fastqcDict['seq_name'] = f"{OrgBID}_{RunID}"
                fastqcDict['orgbatch_id'] = OrgBID
                fastqcDict['run_id'] = RunID
                retLst.append(fdict)
    else:
        fastqcDict = {}
        fastqcDict['seq_name'] = f"{OrgBID}_{RunID}"
        fastqcDict['orgbatch_id'] = OrgBID
        fastqcDict['run_id'] = RunID
        fastqcDict['Seq'] = 'NOT FOUND'
        retLst.append(fdict)
    return(retLst)

# #-----------------------------------------------------------------------------
# def parse_FastQC(AssemblyFolder,OrgBID,RunID):
# #-----------------------------------------------------------------------------
#     BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
#     FastqcDir = os.path.join(BatchDir,"trim","paired","fastqc")
#     outLst = []

#     P1Zip = f"{OrgBID}_{RunID}_P1.zip"
#     with zipfile.ZipFile(os.path.join(FastqcDir,P1Zip)) as z:
#         for filename in z.namelist():
#             if "summary.txt" in filename:
#                 pDict ={}
#                 pDict['Seq'] = 'P1'
#                 with z.open(filename,'r') as f:
#                     reader = csv.reader(io.TextIOWrapper(f),delimiter='\t')
#                     for row in reader:
#                         pDict[row[1]] = row[0]
#                     outLst.append(pDict)

#     P2Zip = f"{OrgBID}_{RunID}_P2.zip"
#     with zipfile.ZipFile(os.path.join(FastqcDir,P2Zip)) as z:
#         for filename in z.namelist():
#             if "summary.txt" in filename:
#                 pDict ={}
#                 pDict['Seq'] = 'P2'
#                 with z.open(filename,'r') as f:
#                     reader = csv.reader(io.TextIOWrapper(f),delimiter='\t')
#                     for row in reader:
#                         pDict[row[1]] = row[0]
#                     outLst.append(pDict)

#     return(outLst)


#-----------------------------------------------------------------------------
def get_CheckM_Info(AssemblyFolder,OrgBID,RunID, Assemblies = ['spades','shovill'], outType = 'contigs', Contamination_cutOff = 5.0):
#-----------------------------------------------------------------------------
    retLst = []

    for Assembly in Assemblies:
        AssInfo = os.path.join(AssemblyFolder,"assembly",Assembly,"checkm",f"info_checkm_{outType}.csv")
        if os.path.exists(AssInfo):

            with open(AssInfo) as file:
                csv_file = csv.DictReader(file)
                binDict = {}
                for line in csv_file:
                    if 'type' in line:
                        if line['type'] == outType:
                            checkmcDict = line
                            checkmcDict['seq_name'] = f"{OrgBID}_{RunID}"
                            checkmcDict['orgbatch_id'] = OrgBID
                            checkmcDict['run_id'] = RunID
                            checkmcDict['assembly'] = Assembly
                            checkmcDict["assembly_qc"] = line['status']
                            # if float(checkmcDict["contamination"]) >= Contamination_cutOff:
                            #     checkmcDict["assembly_qc"] = 'Failed'
                            # else:
                            #     checkmcDict["assembly_qc"] = 'Pass'
                            retLst.append(checkmcDict)
    return(retLst)

# #-----------------------------------------------------------------------------
# def parse_CheckM(AssemblyFolder,OrgBID,RunID,outType="scaffolds"):
# #-----------------------------------------------------------------------------
#     BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
#     CheckmDir = os.path.join(BatchDir,"assembly","spades","checkM","lineage","storage")
#     outDict = {}
#     if os.path.exists(CheckmDir):
#         with open(os.path.join(CheckmDir,"bin_stats_ext.tsv")) as file:
#             tsv_file = csv.reader(file,delimiter="\t")
#             binDict = {}
#             for line in tsv_file:
#                 binDict[line[0]]=eval(line[1])
        
#             outDict["marker_lineage"] = binDict[outType]["marker lineage"]
#             outDict["n_genomes"] = binDict[outType]["# genomes"]
#             outDict["n_predit_genes"] = binDict[outType]["# predicted genes"]
#             outDict["n_markers"] = binDict[outType]["# markers"]
#             outDict["n_marker_sets"] = binDict[outType]["# marker sets"]
#             outDict["n_scaffolds"] = binDict[outType]["# scaffolds"]
#             outDict["genome_size"] = binDict[outType]["Genome size"]
#             outDict["coding_density"] = binDict[outType]["Coding density"]
#             outDict["completeness"] = binDict[outType]["Completeness"]
#             outDict["contamination"] = binDict[outType]["Contamination"]
#             outDict["genome_size"] = binDict[outType]["Genome size"]
#             outDict["gc"] = binDict[outType]["GC"]
#             outDict["gc_std"] = binDict[outType]["GC std"]
#             outDict["n_ambig_bases"] = binDict[outType]["# ambiguous bases"]
#             outDict["longest_scaffold"] = binDict[outType]["Longest scaffold"]
#             outDict["mean_scaffols"] = binDict[outType]["Mean scaffold length"]
#             outDict["N50"] = binDict[outType]["N50 (scaffolds)"]
#             outDict["trans_table"] = binDict[outType]["Translation table"]
#     return(outDict)


#-----------------------------------------------------------------------------
def get_Kraken_Info(FastAFolder,OrgBID,RunID, inType = 'fasta',outType="S",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    KrakenDir = os.path.join(FastAFolder,"kraken")
    KrakenF = os.path.join(KrakenDir,f"{OrgBID}_{inType}.report")
    outLst = []
    if os.path.exists(KrakenF):
        with open(KrakenF) as file:
            tsv_file = csv.reader(file,delimiter="\t")
            for line in tsv_file:
                if line[3] == outType:
                    pctSeq = float(line[0])
                    if pctSeq >= pctCutOff:
                        org_name = line[5].strip()
                        outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
    return(outLst)

# #-----------------------------------------------------------------------------
# def parse_Kraken(AssemblyFolder,OrgBID,RunID,outType="S",inType="fasta",pctCutOff=1.0):
# #-----------------------------------------------------------------------------
#     BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
#     KrakenDir = os.path.join(BatchDir,"kraken")
#     KrakenF = f"{OrgBID}_{RunID}_{inType}.report"
#     outLst = []
#     if os.path.exists(KrakenDir):
#         with open(os.path.join(KrakenDir,KrakenF)) as file:
#             tsv_file = csv.reader(file,delimiter="\t")
#             for line in tsv_file:
#                 if line[3] == outType:
#                     pctSeq = float(line[0])
#                     if pctSeq >= pctCutOff:
#                         org_name = line[5].strip()
#                         outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
#     return(outLst)

# #-----------------------------------------------------------------------------
# def parse_Abricate(AssemblyFolder,OrgBID,RunID,inType="fasta",pctCutOff=1.0):
# #-----------------------------------------------------------------------------
#     BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
#     AbricateDir = os.path.join(BatchDir,"abricate")
#     AbricateF = f"{OrgBID}_{RunID}_{inType}.out"
#     outLst = []
#     if os.path.exists(AbricateDir):
#         with open(os.path.join(AbricateDir,AbricateF)) as file:
#             tsv_file = csv.reader(file,delimiter="\t")
#             for line in tsv_file:
#                 if line[3] == outType:
#                     pctSeq = float(line[0])
#                     if pctSeq >= pctCutOff:
#                         org_name = line[5].strip()
#                         outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
#     return(outLst)

#-----------------------------------------------------------------------------
def get_MLST_Info(FastAFolder,OrgBID,RunID,inType="fasta",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    MlstDir = os.path.join(FastAFolder,"mlst")
    MlstF = os.path.join(MlstDir,f"{OrgBID}_{inType}_mlst.tsv")
    outLst = []
    if os.path.exists(MlstF):
        with open(MlstF) as file:
            tsv_file = csv.reader(file,delimiter="\t",)
            for line in tsv_file:
                outLst.append({'mlst_scheme': line[1], 'mlst_seqtype': line[2], 'mlst_alleles': ";".join(line[3:])})
    return(outLst)

#-----------------------------------------------------------------------------
def get_GTDBTK_Info(FastAFolder,OrgBID,RunID,inType="fasta",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    GTDir = os.path.join(FastAFolder,"gtdbtk")
    GTF = os.path.join(GTDir,f"{OrgBID}_{inType}_gtdbtk.tsv")
    outLst = []
    if os.path.exists(GTF):
        with open(GTF) as file:
            tsv_file = csv.reader(file,delimiter="\t",)
            for line in tsv_file:
                if 'user_genome' not in line:
                    xClass = line[1].split(";")[-1]
                    outLst.append({'gtdbtk_class': xClass, 
                                'gtdbtk_fastani_ref': line[2], 
                                'gtdbtk_fastani_rad': line[3],
                                'gtdbtk_fastani_ani': line[5],
                                'gtdbtk_fastani_af': line[6],
                                })
    return(outLst)

#-----------------------------------------------------------------------------
def get_AMRFinder_Info(FastAFolder,OrgBID,RunID,inType="fasta",pctCutOff=1.0):
#-----------------------------------------------------------------------------
    AmrFinderDir = os.path.join(FastAFolder,"amrfinder")
    AmrFinderF = os.path.join(AmrFinderDir,f"{OrgBID}_{inType}_amrfinder.tsv")
    outLst = []
    if os.path.exists(AmrFinderF):
        with open(AmrFinderF) as file:
            tsv_file = csv.reader(file,delimiter="\t",)
            for line in tsv_file:
                if 'Start' not in line:
                    outLst.append({'amrfinder_contigid': line[1], 
                                'gene_code': line[5], 
                                'gene_name': line[6], 
                                'gene_type': line[8], 
                                'gene_subtype': line[9], 
                                'amr_class': line[10], 
                                'amr_subclass': line[11], 
                                'amrfinder_coverage': line[15],
                                'amrfinder_identitiy': line[16],
                                'amrfinder_closest': line[18],
                                'amrfinder_closestname': line[19],
                                'amr_method':'AMR Finder'
                                })
    return(outLst)

# #-----------------------------------------------------------------------------
# def parse_MLST(AssemblyFolder,OrgBID,RunID,outType="S",inType="fastq",pctCutOff=1.0):
# #-----------------------------------------------------------------------------
#     BatchDir = os.path.join(AssemblyFolder,OrgBID,RunID)
#     KrakenDir = os.path.join(BatchDir,"mlst")
#     KrakenF = f"{OrgBID}_{RunID}_{inType}.report"
#     outLst = []
#     if os.path.exists(KrakenDir):
#         with open(os.path.join(KrakenDir,KrakenF)) as file:
#             tsv_file = csv.reader(file,delimiter="\t")
#             for line in tsv_file:
#                 if line[3] == outType:
#                     pctSeq = float(line[0])
#                     if pctSeq >= pctCutOff:
#                         org_name = line[5].strip()
#                         outLst.append({'org_name': org_name, 'tax_id': line[4], 'pct': pctSeq})
#     return(outLst)




# #-----------------------------------------------------------------------------
# def parse_WGS_COADD(zAssemblyBase):
# #-----------------------------------------------------------------------------

#     lstFastQC = []
#     lstCheckM = []
#     lstKraken = []

#     for dirHead in listFolders(zAssemblyBase):
#         zAssemblyFolder = os.path.join(zAssemblyBase,dirHead)
#         for BatchID in listFolders(zAssemblyFolder):            
#             for RunID in listFolders(os.path.join(zAssemblyFolder,BatchID)):

#                 dirAss = os.path.join(zAssemblyFolder,BatchID,RunID)
#                 print(f"[WGS-COADD] {dirHead} {BatchID} {RunID} --> {dirAss}")

#                 # FastQC
#                 lFastQC=parse_FastQC(zAssemblyFolder,BatchID,RunID)
#                 sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID}
#                 for row in lFastQC:
#                     lstFastQC.append(dict(sDict,**row))

#                 # CheckM
#                 dCheckM = parse_CheckM(zAssemblyFolder,BatchID,RunID)
#                 sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID}
#                 lstCheckM.append(dict(sDict,**dCheckM))

#                 # Kraken
#                 lKraken = parse_Kraken(zAssemblyFolder,BatchID,RunID,outType="S",inType="fasta")
#                 sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}
                
#                 # Abricate
#                 lAbricate = parse_abricate(zAssemblyFolder,BatchID,RunID,outType="S",inType="fasta")
#                 sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}

#                 # MLST
#                 lMlst = parse_Kraken(zAssemblyFolder,BatchID,RunID,outType="S",inType="fastq")
#                 sDict = {'OrgBatch_ID':BatchID, 'Run_ID':RunID, 'ID_Type': 'WGS','ID_Method': 'Kraken2 FastQ','Source':'CO-ADD'}

#                 idLst = []
#                 for i in range(len(lKraken)):
#                     v = lKraken[i]
#                     idLst.append(f"{v['org_name']} ({v['pct']:.1f} pct)")
#                 sDict['ID_Organisms'] = idLst
#                 lstKraken.append(sDict)

#         #         break
#         #     break
#         # break

#     return(lstFastQC,lstCheckM,lstKraken)