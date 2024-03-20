#
import os,sys
import shutil
import pathlib
#
#


#-----------------------------------------------------------------------------
def listFolders(Path):
    return([ name for name in os.listdir(Path) if os.path.isdir(os.path.join(Path, name)) ])
#-----------------------------------------------------------------------------
def listFiles(Path):
    return([ name for name in os.listdir(Path) if os.path.isfile(os.path.join(Path, name)) ])
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def get_RDM(MicroOrgDB):
#-----------------------------------------------------------------------------
    RDM = {
    'rdm' : MicroOrgDB,
    'base' : os.path.join(MicroOrgDB,"Sequence","WGS"),
    'fastq':'01_FastQ',
    'fastq_trim':'01_FastQ_Trim',
    'assembly':'02_Assembly',
    'fasta':'03_FastA',
    'work': 'XX_Work',
    }
    return(RDM)

#-----------------------------------------------------------------------------
def split_BatchID_RunID(batch_run_id):
#-----------------------------------------------------------------------------
    arrStr = batch_run_id.split("_")
    batchID = '_'.join(arrStr[0:3])
    runID = '_'.join(arrStr[3:])
    return batchID, runID

#-----------------------------------------------------------------------------
def get_subdir(OrgBatchID,binsize=200):
#-----------------------------------------------------------------------------
    """
     Gets the SubFolder name based on the XX_NNNN with splits into 200
      GN_0000, GN_0200, GN_0400, ... ,GN_1200, GN_1400, GN_1600
    """
    _org = OrgBatchID.split('_')
    return(f"{_org[0]}_{int(int(_org[1])/binsize)*binsize:04d}")

#-----------------------------------------------------------------------------------
def gen_SeqDict(SeqName,SeqType='WGS'):
#-----------------------------------------------------------------------------------

    OrgBatchID, RunID = split_BatchID_RunID(SeqName)
    OrgID = "_".join(OrgBatchID.split('_')[:2])
    SeqDict = {
        'seq_name'   : f"{OrgBatchID}_{RunID}",
        'orgbatch_id': OrgBatchID,
        'organism_id': OrgID,
        'run_id'     : RunID,
        'sub_dir'    : get_subdir(OrgBatchID),
        'seq_type'   : SeqType,
        'source_code': f"{OrgBatchID}_{RunID}",
        'source_link': f"RDM {SeqType}: {OrgBatchID}_{RunID}",
        'seq_file'   : 'Contigs',
    }
    return(SeqDict)



#-----------------------------------------------------------------------------------
def change_file(rdmDir,oldFile,OrgBatchID,NewBatchID,content=False):
#-----------------------------------------------------------------------------------
    newFile = oldFile.replace(OrgBatchID,NewBatchID)

    print(f"{oldFile} -> {newFile}")
    if os.path.isfile(os.path.join(rdmDir,oldFile)):
        if content:
            print(f" * content {OrgBatchID} -> {NewBatchID}")
            with open(os.path.join(rdmDir,oldFile), 'r') as file:
                filedata = file.read()
            filedata = filedata.replace(OrgBatchID,NewBatchID)
            filedata = filedata.replace(get_subdir(OrgBatchID),get_subdir(NewBatchID))
            with open(os.path.join(rdmDir,newFile), 'w') as file:
                file.write(filedata)
            if oldFile != newFile:
                print(f' * remove {oldFile}')
                os.remove(os.path.join(rdmDir,oldFile))

        elif oldFile != newFile:
            print(f' * move {oldFile} ->{newFile}')
            shutil.move(os.path.join(rdmDir,oldFile), os.path.join(rdmDir,newFile))


#-----------------------------------------------------------------------------------
def rename_SeqName(oldSeqName,newSeqName,MicroOrgDB):
#-----------------------------------------------------------------------------------

    
    RDM = get_RDM(MicroOrgDB)

    oldSeq = gen_SeqDict(oldSeqName)
    newSeq = gen_SeqDict(newSeqName)

    # 01_FastQ 
    # ----------------------------------------------------------
    #     rename file
    # 01_FastQ\oldDir\oldName_R1.fastq.gz 
    # 01_FastQ\oldDir\oldName_R2.fastq.gz 

    fqDir = os.path.join(RDM['base'],RDM['fastq'],oldSeq['sub_dir'])
    if os.path.isdir(fqDir):
        print(RDM['fastq'])
        for fqFiles in listFiles(fqDir):
            if oldSeqName in fqFiles:
                change_file(fqDir,fqFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                print(" - ", fqFiles)
    
    # 01_FastQ_trim
    # ----------------------------------------------------------
    # rename file
    # 01_FastQ_Trim\oldDir\oldName\oldName.log.gz
    # 01_FastQ_Trim\oldDir\oldName\oldName_P1.fastq.gz
    # 01_FastQ_Trim\oldDir\oldName\oldName.P2.fastq.gz
    # 01_FastQ_Trim\oldDir\oldName\fastqc\oldName_P1.fastq.html
    # 01_FastQ_Trim\oldDir\oldName\fastqc\oldName_P1.fastq.zip
    # 01_FastQ_Trim\oldDir\oldName\fastqc\oldName_P2.fastq.html
    # 01_FastQ_Trim\oldDir\oldName\fastqc\oldName_P2.fastq.zip

    fqtDir = os.path.join(RDM['base'],RDM['fastq_trim'],oldSeq['sub_dir'],oldSeq['seq_name'])
    new_fqtDir = os.path.join(RDM['base'],RDM['fastq_trim'],newSeq['sub_dir'],newSeq['seq_name'])

    if os.path.isdir(fqtDir):
        print(RDM['fastq_trim'])
        for fqtFiles in listFiles(fqtDir):
            if oldSeqName in fqtFiles:
                change_file(fqtDir,fqtFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                print(" - ", fqtFiles)

        fqcDir = os.path.join(RDM['base'],RDM['fastq_trim'],oldSeq['sub_dir'],oldSeq['seq_name'],'fastqc')
        for fqcFiles in listFiles(fqcDir):
            if oldSeqName in fqcFiles:
                change_file(fqcDir,fqcFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                print(" - FastQC - ", fqcFiles)
        os.rename(fqtDir, new_fqtDir)

    # 02_Assembly
    # ----------------------------------------------------------
    # rename file 
    # 02_Assembly\oldDir\oldName\shovill\oldName_contigs.fa.gz
    # 02_Assembly\oldDir\oldName\shovill\oldName_contigs.gfa.gz
    # 02_Assembly\oldDir\oldName\shovill\oldName_shovill.correction
    # 02_Assembly\oldDir\oldName\shovill\oldName_shovill_uncorrected_spades.fa.gz
    # 02_Assembly\oldDir\oldName\shovill\checkm\oldName_tab_complecont_stats.txt
    # 02_Assembly\oldDir\oldName\spades\oldName_contigs.fa.gz
    # 02_Assembly\oldDir\oldName\spades\oldName_scaffolds.fa.gz
    # 02_Assembly\oldDir\oldName\spades\checkm\oldName_tab_complecont_stats.txt

    # replace content (oldDir\oldName with newDir\newName):
    # 02_Assembly\oldDir\oldName\info_assembly.json


    assDir = os.path.join(RDM['base'],RDM['assembly'],oldSeq['sub_dir'],oldSeq['seq_name'])
    new_assDir = os.path.join(RDM['base'],RDM['assembly'],newSeq['sub_dir'],newSeq['seq_name'])

    if os.path.isdir(assDir):
        print(RDM['assembly'])
        change_file(assDir,'info_assembly.json',oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=True)

        for assType in listFolders(assDir):
            assTypeDir = os.path.join(assDir,assType)
            assCheckMDir = os.path.join(assTypeDir,'checkm')

            for assFiles in listFiles(assTypeDir):
                if oldSeqName in assFiles:
                    change_file(assTypeDir,assFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                    print(" - ", assFiles)

            for chmFiles in listFiles(assCheckMDir):
                if oldSeqName in chmFiles:
                    change_file(assCheckMDir,chmFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                    print(" - CheckM - ", chmFiles)
        os.rename(assDir, new_assDir)

    # 03_FastA
    # ----------------------------------------------------------
    # rename file
    # 03_FastA\oldDir\oldName\oldName_contigs_filtered.fasta
    # 03_FastA\oldDir\oldName\oldName_<shovill/spade>_checkm_stats.txt

    # rename file and replace content
    # 03_FastA\oldDir\oldName\abricate
    # 03_FastA\oldDir\oldName\amrfinder
    # 03_FastA\oldDir\oldName\gtdbtk
    # 03_FastA\oldDir\oldName\kraken
    # 03_FastA\oldDir\oldName\mlst    

    faDir = os.path.join(RDM['base'],RDM['fasta'],oldSeq['sub_dir'],oldSeq['seq_name'])
    new_faDir = os.path.join(RDM['base'],RDM['fasta'],newSeq['sub_dir'],newSeq['seq_name'])

    if os.path.isdir(faDir):
        print(RDM['fasta'])

        for faType in listFolders(faDir):
            faTypeDir = os.path.join(faDir,faType)

            for amrFiles in listFiles(faTypeDir):
                # only *.tsv or *.report 
                if oldSeq['orgbatch_id'] in amrFiles:
                    chg_content = amrFiles.endswith('.tsv') or amrFiles.endswith('.report')
                    change_file(faTypeDir,amrFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=chg_content)
                    print(" - ", amrFiles,chg_content)

        for faFiles in listFiles(faDir):
            if oldSeqName in faFiles:
                change_file(faDir,faFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)
                print(" - ", faFiles)

        os.rename(faDir, new_faDir)

    # Test
    # ----------------------------------------------------------
    #change_file(assDir,'info_assembly_test.json',oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=True)
    # tfile = "GN_0696_00_fasta_mlst_test.tsv"
    # tDir = os.path.join(RDM['base'],RDM['fasta'],oldSeq['sub_dir'],oldSeq['seq_name'],'mlst')
    # change_file(tDir,tfile,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=False)

# #-----------------------------------------------------------------------------------
def rename_Work(oldSeqName,newSeqName,MicroOrgDB):
# #-----------------------------------------------------------------------------------

    RDM = get_RDM(MicroOrgDB)

    oldSeq = gen_SeqDict(oldSeqName)
    newSeq = gen_SeqDict(newSeqName)

    wDir = os.path.join(RDM['base'],RDM['work'],oldSeq['seq_name'])
    new_wDir = os.path.join(RDM['base'],RDM['work'],newSeq['seq_name'])

    if os.path.isdir(wDir):
        print(RDM['work'])
        for wFiles in listFiles(wDir):
            if oldSeqName in wFiles:
                change_file(wDir,wFiles,oldSeq['orgbatch_id'],newSeq['orgbatch_id'],content=True)
                print(" - ", wFiles)
        os.rename(wDir, new_wDir)

# -----------------------------------------------------------------
def main():
# -----------------------------------------------------------------
    MicroOrgDB = 'I:\MICROORGDB-Q5308'
    # rename_SeqName('GN_0696_00_AGRF_R001','GP_0313_00_AGRF_R001',MicroOrgDB)
    # rename_SeqName('GN_0697_01_AGRF_R001','GP_0314_01_AGRF_R001',MicroOrgDB)
    # rename_SeqName('GN_1110_02_AGRF_R003','GP_0331_02_AGRF_R003',MicroOrgDB)
    # rename_SeqName('GN_1113_01_AGRF_R004','GP_0333_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1178_01_AGRF_R004','GP_0339_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1192_01_AGRF_R004','GP_0340_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1115_03_AGRF_R004','GP_0341_03_AGRF_R004',MicroOrgDB)

    # rename_SeqName('GN_0952_02_AGRF_R004','GN_1382_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_0981_01_AGRF_R002','GN_1383_01_AGRF_R002',MicroOrgDB)
    # rename_SeqName('GN_0982_02_AGRF_R004','GN_1384_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1149_02_AGRF_R003','GN_1385_02_AGRF_R003',MicroOrgDB)
    # rename_SeqName('GN_1254_02_AGRF_R004','GN_1386_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1258_01_AGRF_R004','GN_1387_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1270_03_AGRF_R002','GN_1388_03_AGRF_R002',MicroOrgDB)
    # rename_SeqName('GN_1274_02_AGRF_R002','GN_1389_02_AGRF_R002',MicroOrgDB)
    # rename_SeqName('GN_1275_02_AGRF_R004','GN_1390_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1282_02_AGRF_R004','GN_1391_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1285_02_AGRF_R002','GN_1392_02_AGRF_R002',MicroOrgDB)
    # rename_SeqName('GN_1308_01_AGRF_R004','GN_1393_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1307_01_AGRF_R004','GN_1394_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1316_02_AGRF_R004','GN_1395_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1318_01_AGRF_R004','GN_1396_01_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1323_02_AGRF_R004','GN_1397_02_AGRF_R004',MicroOrgDB)
    # rename_SeqName('GN_1334_02_AGRF_R004','GN_1398_02_AGRF_R004',MicroOrgDB)

    # rename_Work('GN_0696_00_AGRF_R001','GP_0313_00_AGRF_R001',MicroOrgDB)
    # rename_Work('GN_0697_01_AGRF_R001','GP_0314_01_AGRF_R001',MicroOrgDB)
    # rename_Work('GN_1110_02_AGRF_R003','GP_0331_02_AGRF_R003',MicroOrgDB)
    # rename_Work('GN_1113_01_AGRF_R004','GP_0333_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1178_01_AGRF_R004','GP_0339_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1192_01_AGRF_R004','GP_0340_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1115_03_AGRF_R004','GP_0341_03_AGRF_R004',MicroOrgDB)

    # rename_Work('GN_0952_02_AGRF_R004','GN_1382_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_0981_01_AGRF_R002','GN_1383_01_AGRF_R002',MicroOrgDB)
    # rename_Work('GN_0982_02_AGRF_R004','GN_1384_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1149_02_AGRF_R003','GN_1385_02_AGRF_R003',MicroOrgDB)
    # rename_Work('GN_1254_02_AGRF_R004','GN_1386_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1258_01_AGRF_R004','GN_1387_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1270_03_AGRF_R002','GN_1388_03_AGRF_R002',MicroOrgDB)
    # rename_Work('GN_1274_02_AGRF_R002','GN_1389_02_AGRF_R002',MicroOrgDB)
    # rename_Work('GN_1275_02_AGRF_R004','GN_1390_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1282_02_AGRF_R004','GN_1391_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1285_02_AGRF_R002','GN_1392_02_AGRF_R002',MicroOrgDB)
    # rename_Work('GN_1308_01_AGRF_R004','GN_1393_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1307_01_AGRF_R004','GN_1394_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1316_02_AGRF_R004','GN_1395_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1318_01_AGRF_R004','GN_1396_01_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1323_02_AGRF_R004','GN_1397_02_AGRF_R004',MicroOrgDB)
    # rename_Work('GN_1334_02_AGRF_R004','GN_1398_02_AGRF_R004',MicroOrgDB)


#==============================================================================
if __name__ == "__main__":
    print("-------------------------------------------------------------------")
    print("Running : ",sys.argv)
    print("-------------------------------------------------------------------")
    main()
    print("...................................................................")
#==============================================================================
